// ============================================================
// Nuclei Binary Mask Generator
// Fiji/ImageJ Macro
//
// Opens one or more nucleus images selected by the user,
// applies automatic thresholding to produce a binary mask,
// and saves each mask as a TIFF alongside a log of settings.
//
// Workflow per image:
//   1. Convert to 8-bit greyscale
//   2. (Optional) Gaussian blur to reduce noise
//   3. Auto-threshold using the chosen algorithm
//   4. Convert to binary mask (0 / 255)
//   5. (Optional) Fill holes and/or remove small debris
//   6. Save mask as <originalName>_mask.tif
// ============================================================

// ---------------------------------------------------------------
// Configuration — edit these to suit your images
// ---------------------------------------------------------------

// Auto-threshold algorithm.
// Common choices for DAPI/Hoechst nuclei:
//   "Default"  "Otsu"  "Triangle"  "Li"  "Mean"  "RenyiEntropy"
var THRESHOLD_METHOD = "Otsu";

// Apply a Gaussian blur before thresholding to reduce noise?
var APPLY_BLUR = true;
var BLUR_SIGMA = 2.0;   // pixels; increase for noisier images

// Fill holes inside detected nuclei?
var FILL_HOLES = true;

// Remove small objects (salt noise / debris)?
var REMOVE_SMALL = true;
var MIN_SIZE_PX  = 50;  // objects smaller than this (pixels²) are removed

// Suffix appended to the original filename for the mask output
var MASK_SUFFIX = "_mask";

// ---------------------------------------------------------------
// Entry point
// ---------------------------------------------------------------
main();

function main() {
	// Let the user pick one or more image files
	// (hold Shift or Ctrl to multi-select in the dialog)
	imagePaths = File.openDialog("Select a nucleus image (hold Shift/Ctrl for multiple)");

	// File.openDialog returns a single path as a string.
	// For multi-file selection we need the newer OpenDialog approach.
	// We therefore offer a folder-mode fallback if the user wants batch.
	choice = getBoolean("Process a single file?\n \n" +
	                    "Click YES to process ONE image you just chose.\n" +
	                    "Click NO  to process ALL images in a FOLDER.");

	if (choice) {
		// ---- Single file mode ----
		if (imagePaths == "") exit("No file selected. Macro aborted.");
		outputDir = File.getParent(imagePaths) + File.separator;
		processImage(imagePaths, outputDir);
	} else {
		// ---- Folder (batch) mode ----
		inputDir = getDirectory("Select the folder containing nucleus images");
		if (inputDir == "") exit("No folder selected. Macro aborted.");

		outputDir = getDirectory("Select the OUTPUT folder for mask TIFFs");
		if (outputDir == "") exit("No output folder selected. Macro aborted.");

		batchProcess(inputDir, outputDir);
	}

	showMessage("Done", "Mask generation complete.\nCheck the Log window for details.");
}

// ---------------------------------------------------------------
// Batch: process every image in a folder (non-recursive)
// Supported extensions: tif, tiff, png, jpg, jpeg, czi, lif, nd2
// ---------------------------------------------------------------
function batchProcess(inputDir, outputDir) {
	extensions = newArray(".tif", ".tiff", ".png", ".jpg", ".jpeg",
	                      ".czi", ".lif", ".nd2", ".bmp");

	fileList = getFileList(inputDir);
	found    = 0;

	for (i = 0; i < fileList.length; i++) {
		name = fileList[i];
		if (File.isDirectory(inputDir + name)) continue;

		lower = toLowerCase(name);
		isImage = false;
		for (e = 0; e < extensions.length; e++) {
			if (endsWith(lower, extensions[e])) { isImage = true; break; }
		}

		if (isImage) {
			processImage(inputDir + name, outputDir);
			found++;
		}
	}

	if (found == 0) print("WARNING: No recognised image files found in " + inputDir);
}

// ---------------------------------------------------------------
// Core: open one image, threshold, save mask
// ---------------------------------------------------------------
function processImage(filePath, outputDir) {
	fileName = File.getName(filePath);
	print("Processing: " + filePath);

	// Open via Bio-Formats so exotic formats (CZI etc.) also work
	run("Bio-Formats Importer",
	    "open=[" + filePath + "] " +
	    "autoscale color_mode=Grayscale view=Hyperstack stack_order=XYCZT " +
	    "series_list=1");

	if (nImages == 0) {
		print("  ERROR: Could not open " + filePath);
		return;
	}

	imgID = getImageID();
	selectImage(imgID);

	// ---- 1. Flatten to single 2-D greyscale plane ----
	// If the image is a stack/hyperstack, max-project across Z
	getDimensions(width, height, channels, slices, frames);

	if (slices > 1) {
		run("Z Project...", "projection=[Max Intensity]");
		close("\\Others");   // close the original stack
		imgID = getImageID();
		selectImage(imgID);
		print("  Z-projected " + slices + " slices (Max Intensity)");
	}

	// If multi-channel, ask the user which channel to use
	getDimensions(width, height, channels, slices, frames);
	if (channels > 1) {
		channelNum = getNumber("Image has " + channels + " channels.\n" +
		                       "Which channel contains the nuclei? (1-" + channels + ")",
		                       1);
		Stack.setChannel(channelNum);
		run("Duplicate...", "title=nuclei_channel duplicate channels=" + channelNum);
		close("\\Others");
		imgID = getImageID();
		selectImage(imgID);
		print("  Extracted channel " + channelNum);
	}

	// Convert to 8-bit
	run("8-bit");

	// ---- 2. Gaussian blur (noise suppression) ----
	if (APPLY_BLUR) {
		run("Gaussian Blur...", "sigma=" + BLUR_SIGMA);
		print("  Gaussian blur sigma=" + BLUR_SIGMA);
	}

	// ---- 3. Auto-threshold ----
	setAutoThreshold(THRESHOLD_METHOD + " dark");
	run("Convert to Mask");
	print("  Threshold method: " + THRESHOLD_METHOD);

	// ---- 4. Fill holes ----
	if (FILL_HOLES) {
		run("Fill Holes");
		print("  Holes filled");
	}

	// ---- 5. Remove small objects ----
	if (REMOVE_SMALL) {
		// Analyse particles: exclude objects smaller than MIN_SIZE_PX,
		// redirect to nothing, output a clean binary mask
		run("Analyze Particles...",
		    "size=" + MIN_SIZE_PX + "-Infinity " +
		    "circularity=0.00-1.00 show=Masks exclude clear");
		// "Mask of ..." window is now the foreground; close the raw mask
		maskTitle = getTitle();   // "Mask of <originalTitle>"
		// Bring the particle mask to the front (it becomes the active image)
		imgID = getImageID();
		selectImage(imgID);
		// Make it a proper 8-bit binary (0/255)
		setMinAndMax(0, 255);
		run("8-bit");
		print("  Small objects < " + MIN_SIZE_PX + " px² removed");
	}

	// ---- 6. Build output path and save ----
	// Strip known extensions from the filename
	baseName = fileName;
	knownExts = newArray(".tif",".tiff",".png",".jpg",".jpeg",".czi",".lif",".nd2",".bmp");
	for (e = 0; e < knownExts.length; e++) {
		if (endsWith(toLowerCase(baseName), knownExts[e])) {
			baseName = substring(baseName, 0, lengthOf(baseName) - lengthOf(knownExts[e]));
			break;
		}
	}

	outName = baseName + MASK_SUFFIX + ".tif";
	outPath = outputDir + outName;

	saveAs("Tiff", outPath);
	print("  Saved mask: " + outPath);

	close();
}
