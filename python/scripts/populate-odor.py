#!/usr/local/bin/python3
from pipeline import experiment, reso, meso, fuse, stack, pupil, treadmill

# # Scans
# for priority in range(120, -130, -10):  # highest to lowest priority
#     next_scans = (experiment.AutoProcessing() & 'priority > {}'.format(priority) &
#                   (experiment.Scan() & 'scan_ts > "2019-01-01 00:00:00"'))

next_scans = (experiment.AutoProcessing  &
              (experiment.Scan & 'scan_ts > "2019-01-01 00:00:00"'))

kwargs = dict(suppress_errors=True, display_progress=True)


# treadmill, pupil
treadmill.Treadmill.populate(next_scans, **kwargs)
pupil.Eye.populate(next_scans, **kwargs)
pupil.FittedPupil.populate(next_scans, **kwargs)

# stack
stack.StackInfo.populate(stack.CorrectionChannel, **kwargs)
stack.Quality.populate(**kwargs)
stack.RasterCorrection.populate(**kwargs)
stack.MotionCorrection.populate(**kwargs)
stack.Stitching.populate(**kwargs)
stack.CorrectedStack.populate(**kwargs)

# reso/meso
for pipe in [reso, meso]:
    pipe.ScanInfo.populate(next_scans, **kwargs)
    pipe.Quality.populate(next_scans, **kwargs)
    pipe.RasterCorrection.populate(next_scans, **kwargs)
    pipe.MotionCorrection.populate(next_scans, **kwargs)
    pipe.SummaryImages.populate(next_scans, **kwargs)
    pipe.Segmentation.populate(next_scans, **kwargs)
    pipe.Fluorescence.populate(next_scans, **kwargs)
    pipe.MaskClassification.populate(next_scans, {'classification_method': 2},
                                     **kwargs)
    pipe.ScanSet.populate(next_scans, **kwargs)
    pipe.Activity.populate(next_scans, {'spike_method': 5}, reserve_jobs=True,
                           suppress_errors=True)
    full_scans = (pipe.ScanInfo.proj() & pipe.Activity) - (pipe.ScanInfo.Field -
                                                           pipe.Activity)
    pipe.ScanDone.populate(full_scans & next_scans, reserve_jobs=True,
                           suppress_errors=True)

# fuse
fuse.MotionCorrection.populate(next_scans, **kwargs)
fuse.ScanSet.populate(next_scans, **kwargs)
fuse.Activity.populate(next_scans, **kwargs)
fuse.ScanDone.populate(next_scans, **kwargs)

# more stack (needs corrected fields)
stack.PreprocessedStack.populate(stack.SegmentationTask, reserve_jobs=True,
                                 suppress_errors=True)
stack.FieldSegmentation.populate(**kwargs)
stack.PreprocessedStack.populate(stack.RegistrationTask.proj(session='stack_session',
                                                             channel='stack_channel'),
                                 **kwargs)
stack.Registration.populate(**kwargs)

# # tune (these are memory intensive)
# tune.STA.populate(next_scans, **kwargs)
# tune.STAQual.populate(next_scans, **kwargs)
# tune.STAExtent.populate(next_scans, **kwargs)
#
# tune.CaMovie.populate(next_scans, **kwargs)
# tune.Drift.populate(next_scans, **kwargs)
# tune.OriDesign.populate(next_scans, **kwargs)
# tune.OriMap.populate(next_scans, **kwargs)
# tune.Cos2Map.populate(next_scans, **kwargs)
# tune.OriMapQuality.populate(next_scans, **kwargs)
#
# # tune.OracleMap.populate(next_scans, **kwargs)
# # tune.MovieOracle.populate(next_scans, **kwargs)
# # tune.MovieOracleTimeCourse.populate(next_scans, **kwargs)
#
# tune.CaTimes.populate(next_scans, **kwargs)
# tune.Ori.populate(next_scans, **kwargs)
