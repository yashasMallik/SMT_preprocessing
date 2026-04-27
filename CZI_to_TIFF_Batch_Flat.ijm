// ============================================================
// CZI to TIFF Batch Converter (Recursive, Flat Output)
// Fiji/ImageJ Macro
//
// Converts all .czi files in a selected folder (and all
// sub-folders) to .tif files saved in a SINGLE flat output
// folder. Output filenames encode the relative folder path so
// nothing is overwritten and the origin is always clear.
//
// Example:
//   Input:  RootFolder/Condition_A/Day1/sample.czi
//   Output: OutputFolder/Condition_A_Day1_sample.tif
// ============================================================

// --- Configuration ---
// Separator used between folder names and the file name in the output filename.
// "_" gives:  Condition_A_Day1_sample.tif
// "--" gives: Condition_A--Day1--sample.tif
var SEPARATOR = "_";

// Set to true  → saves each scene in a multi-scene CZI as a separate TIFF
// Set to false → saves only the first series/scene
var SPLIT_SCENES = true;

// Set to true  → skip conversion if the output TIFF already exists
var SKIP_EXISTING = true;

// ---------------------------------------------------------------
// Entry point
// ---------------------------------------------------------------
main();

function main() {
	inputDir = getDirectory("Select the ROOT folder containing CZI files");
	if (inputDir == "") exit("No input folder selected. Macro aborted.");

	outputDir = getDirectory("Select the OUTPUT folder for TIFF files");
	if (outputDir == "") exit("No output folder selected. Macro aborted.");

	// Normalise: ensure trailing separator
	if (!endsWith(inputDir,  File.separator)) inputDir  = inputDir  + File.separator;
	if (!endsWith(outputDir, File.separator)) outputDir = outputDir + File.separator;

	convertedCount = 0;
	skippedCount   = 0;
	errorCount     = 0;

	processFolder(inputDir, inputDir, outputDir);

	showMessage("Batch Conversion Complete",
		"Converted : " + convertedCount + " file(s)\n" +
		"Skipped   : " + skippedCount   + " file(s)\n" +
		"Errors    : " + errorCount      + " file(s)");
}

// ---------------------------------------------------------------
// Recursively walk folders; collect CZI files
// ---------------------------------------------------------------
function processFolder(currentDir, inputRootDir, outputDir) {
	fileList = getFileList(currentDir);

	for (i = 0; i < fileList.length; i++) {
		name     = fileList[i];
		fullPath = currentDir + name;

		if (File.isDirectory(fullPath)) {
			processFolder(fullPath + File.separator, inputRootDir, outputDir);

		} else if (endsWith(toLowerCase(name), ".czi")) {
			// Build a flat prefix from the relative path of the parent folder.
			// e.g. "Condition_A/Day1/" → "Condition_A_Day1"
			relativeDir = substring(currentDir, lengthOf(inputRootDir));
			prefix      = buildPrefix(relativeDir);

			convertCZI(fullPath, outputDir, name, prefix);
		}
	}
}

// ---------------------------------------------------------------
// Turn a relative directory path into a flat filename prefix.
// "Condition_A/Day1/"  →  "Condition_A_Day1"  (using SEPARATOR)
// ""                   →  ""  (file is directly in the root)
// ---------------------------------------------------------------
function buildPrefix(relativeDir) {
	// Strip trailing separator
	if (endsWith(relativeDir, File.separator)) {
		relativeDir = substring(relativeDir, 0, lengthOf(relativeDir) - 1);
	}
	if (relativeDir == "") return "";

	// Replace all path separators with our chosen SEPARATOR
	prefix = replace(relativeDir, File.separator, SEPARATOR);
	// Also replace any remaining / or \ (cross-platform safety)
	prefix = replace(prefix, "/",  SEPARATOR);
	prefix = replace(prefix, "\\", SEPARATOR);
	// Collapse runs of the separator that might have appeared
	while (indexOf(prefix, SEPARATOR + SEPARATOR) >= 0) {
		prefix = replace(prefix, SEPARATOR + SEPARATOR, SEPARATOR);
	}
	return prefix;
}

// ---------------------------------------------------------------
// Convert a single CZI file and save to the flat output folder.
// ---------------------------------------------------------------
function convertCZI(cziPath, outputDir, fileName, prefix) {
	// Strip extension (case-insensitive)
	baseName = fileName;
	if (endsWith(toLowerCase(baseName), ".czi")) {
		baseName = substring(baseName, 0, lengthOf(baseName) - 4);
	}

	// Build the flat output base name:  prefix_baseName  or just baseName
	if (prefix != "") {
		flatBase = prefix + SEPARATOR + baseName;
	} else {
		flatBase = baseName;
	}

	// Quick existence check before opening the file
	firstOutput = outputDir + flatBase + ".tif";
	if (SKIP_EXISTING && File.exists(firstOutput)) {
		print("Skipped (already exists): " + cziPath);
		skippedCount++;
		return;
	}

	print("Converting: " + cziPath);
	print("  Output prefix: " + flatBase);

	importOptions = "open=[" + cziPath + "] " +
	                "autoscale " +
	                "color_mode=Default " +
	                "rois_import=[ROI manager] " +
	                "view=Hyperstack " +
	                "stack_order=XYCZT";

	run("Bio-Formats Macro Extensions");
	Ext.setId(cziPath);
	Ext.getSeriesCount(seriesCount);

	if (seriesCount <= 0) {
		print("  ERROR: Could not read series count for " + cziPath);
		errorCount++;
		Ext.close();
		return;
	}

	print("  Series/scenes found: " + seriesCount);

	maxSeries = seriesCount;
	if (!SPLIT_SCENES) maxSeries = 1;

	for (s = 0; s < maxSeries; s++) {
		seriesOptions = importOptions + " series_list=" + (s + 1);
		run("Bio-Formats Importer", seriesOptions);

		if (nImages == 0) {
			print("  WARNING: No image opened for series " + (s + 1));
			errorCount++;
			continue;
		}

		imgID = getImageID();
		selectImage(imgID);

		// Build output filename
		if (seriesCount == 1 || !SPLIT_SCENES) {
			outName = flatBase + ".tif";
		} else {
			Ext.setSeries(s);
			Ext.getSeriesName(seriesName);
			seriesName = replace(seriesName, "[^a-zA-Z0-9_\\-]", "_");
			if (seriesName == "" || seriesName == "_") {
				seriesName = "scene_" + (s + 1);
			}
			outName = flatBase + SEPARATOR + seriesName + ".tif";
		}

		outPath = outputDir + outName;
		saveAs("Tiff", outPath);
		print("  Saved: " + outPath);

		close();
		convertedCount++;
	}

	Ext.close();
}
