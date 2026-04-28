TrackDBET_2<-createTrackll(folder=DBET_2, input=2, cores=1, frameRecord=T)
TrackDBET_6<-createTrackll(folder=DBET_6, input=2, cores=1, frameRecord=T)
TrackDBET_20<-createTrackll(folder=DBET_20, input=2, cores=1, frameRecord=T)
TrackMZ_2<-createTrackll(folder=MZ_2, input=2, cores=1, frameRecord=T)
TrackMZ_6<-createTrackll(folder=MZ_6, input=2, cores=1, frameRecord=T)
TrackMZ_20<-createTrackll(folder=MZ_20 , input=2, cores=1, frameRecord=T)


TrackDBET_2.fi<-filterTrack(trackll=TrackDBET_2, filter=c(min=6,max=Inf))
TrackDBET_6.fi<-filterTrack(trackll=TrackDBET_6, filter=c(min=6,max=Inf))
TrackDBET_20.fi<-filterTrack(trackll=TrackDBET_20, filter=c(min=6,max=Inf))
TrackMZ_2.fi<-filterTrack(trackll=TrackMZ_2, filter=c(min=6,max=Inf))
TrackMZ_6.fi<-filterTrack(trackll=TrackMZ_6, filter=c(min=6,max=Inf))
TrackMZ_20.fi<-filterTrack(trackll=TrackMZ_20, filter=c(min=6,max=Inf))


TrackDBET_2.fi.li<-linkSkippedFrames(trackll = TrackDBET_2.fi,tolerance=2,maxSkip =5,cores=1)
TrackDBET_6.fi.li<-linkSkippedFrames(trackll = TrackDBET_6.fi,tolerance=2,maxSkip =5,cores=1)
TrackDBET_20.fi.li<-linkSkippedFrames(trackll = TrackDBET_20.fi,tolerance=2,maxSkip =5,cores=1)
TrackMZ_2.fi.li<-linkSkippedFrames(trackll = TrackMZ_2.fi,tolerance=2,maxSkip =5,cores=1)
TrackMZ_6.fi.li<-linkSkippedFrames(trackll = TrackMZ_6.fi,tolerance=2,maxSkip =5,cores=1)
TrackMZ_20.fi.li<-linkSkippedFrames(trackll = TrackMZ_20.fi,tolerance=2,maxSkip =5,cores=1)


TrackDBET_2.fi.li.ma <- maskTracks(folder = dbet2, trackll = TrackDBET_2.fi.li)
TrackDBET_6.fi.li.ma <- maskTracks(folder = dbet6, trackll = TrackDBET_6.fi.li)
TrackDBET_20.fi.li.ma <- maskTracks(folder = dbet20, trackll = TrackDBET_20.fi.li)
TrackMZ_2.fi.li.ma <- maskTracks(folder = mz2, trackll = TrackMZ_2.fi.li)
TrackMZ_6.fi.li.ma <- maskTracks(folder = mz6, trackll = TrackMZ_6.fi.li)
TrackMZ_20.fi.li.ma <- maskTracks(folder = mz20, trackll = TrackMZ_20.fi.li)

//not used
TrackDBET_2.fi.li.ma.me <-mergeTracks(dbet2, TrackDBET_2.fi.li.ma)
TrackDBET_6.fi.li.ma.me <-mergeTracks(dbet6, TrackDBET_6.fi.li.ma)
TrackDBET_20.fi.li.ma.me <-mergeTracks(dbet20, TrackDBET_20.fi.li.ma)
TrackMZ_2.fi.li.ma.me <-mergeTracks(mz2, TrackMZ_2.fi.li.ma)
TrackMZ_6.fi.li.ma.me <-mergeTracks(mz6, TrackMZ_6.fi.li.ma)
TrackMZ_20.fi.li.ma.me <-mergeTracks(mz20, TrackMZ_20.fi.li.ma)



fitRT(trackll= TrackDBET_2.fi.li.ma, x.max=100, N.min=1.5, t.interval=0.2)
fitRT(trackll= TrackDBET_6.fi.li.ma, x.max=100, N.min=1.5, t.interval=0.2)
fitRT(trackll= TrackDBET_20.fi.li.ma, x.max=100, N.min=1.5, t.interval=0.2)
fitRT(trackll= TrackMZ_2.fi.li.ma, x.max=100, N.min=1.5, t.interval=0.2)
fitRT(trackll= TrackMZ_6.fi.li.ma, x.max=100, N.min=1.5, t.interval=0.2)
fitRT(trackll= TrackMZ_20.fi.li.ma, x.max=100, N.min=1.5, t.interval=0.2)
