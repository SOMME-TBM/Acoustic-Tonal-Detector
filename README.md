# Acoustic-Tonal-Detector
An automatic tonal detector combining energy detection with a contour shape fitting algorithm for passive acoustic monitoring.

This detector is described in the following paper: Caron Delbosc N., Mathias D., Chauvaud L., Ahonen H., Llobet S.M., Lydersen C., Kovacs K.M.  and  Richard G. (2025) Development of a conservative automated tonal detector with high performance at a large temporal scale, Bioacoustics, DOI: 10.1080/09524622.2025.2568528
Please cite this article the algorithm is used.

**ABSTRACT**
In the context of global warming, monitoring changes in marine fauna distribution is crucial, especially in remote regions. Passive acoustic monitoring (PAM) has become a valuable method to detect vocal species over long periods, particularly in Arctic waters. However, the growing volume of data collected raises new analytical challenges. In this study, we developed an automatic tonal detector combining energy detection with a contour shape fitting algorithm, tested on bearded seal (Erignathus barbatus) vocalisations. This Arctic species, highly vocal during the breeding season and living in remote areas, is well-suited to PAM-based monitoring. We assessed the detector’s performance across areas and seasons. While the detector showed high accuracy in identifying individual vocalisations (100% correct detections), it also missed a substantial number of calls (77.5% false negatives). Despite this conservative detection approach, it effectively captured seasonal presence patterns, with 80% to 96% agreement with manual annotations in Svalbard showing that this low recall rate does not compromise the reliability of the study in an ecological approach. This tool offers a rapid, first-order estimation method for detecting bearded seals - or other tonal sound-producing species - from large acoustic datasets collected over extended time periods.


**Tonal detector - description**

This tonal detector finds energy peaks over the signal and then assessed  tonal contour.

The first step is an energy detector, which estimates the average sound level exposure (aPSD) as the power spectral density level integrated over a frequency bandwidths (bandwidth_call_freq) for a given duration (dt_aPSD). The detector considers the distribution of aPSD values over the signal with the hypothesis that the distribtuion should be normal, and thus based upon the mode, median or mean of aPSD values set a threshold above which values are considered as overdisperes (i.e. few occurence of high energy):

1- First step doesn't discriminate energetic sounds, thus a first selection is based upon the tonal durations with a minimum and a maximum  duration possible for the tonals (window_call_durations)
2- The second step is a contour-based detector. The contour (i.e. frequency) is recovered from the spectrogram of the potential tonals detected based upon the energy criteria, by finding which frequency has the maximum energy for every time unit of the spectrogram. From tonals definition, the contour should be smooth and thus a fit model could be applied (fit_type). Then, a R² threshold is set to considered whether the energetic sound follows accurately a descriptive fit model.


**Detector output**

'raw_tonals' is a cell array with all contours detected and extracted:
every cell (tonal) is a double-precision array [time frequency]
the time is in second from the begininig of the signal (sig_acou)

'smooth_tonals' is the same as raw_tonals but with smoothed contours (using a pchip relationship)

'energetic_detector' is a matrix (type double) with two column:
First col: time in second when the average PSD was found above a threshold (resolution of this column is at 'dt_aPSD' seconds)
Second col: numer of the event (several line with the same number = same event

'energetic_duration_detector' is a matrix (type double) with three column derived from 'energetic_detector' by considering onmy events long enough and not too long (i.e within [window_call_durations]):
1st col: begining time in second of the event kept from 'energetic_detector'
2nd col: end time in second of the event kept from 'energetic_detector'
3rd col: number of the event (kept from 'energetic_detector')



**Detector parameters  - INPUTS**

sig_acou: acoustic signal
fs: sampling rate of the acoustic signal

dt_aPSD : aPSD resoloution => time window over which thre PSD is estimated to then average between the bandwidth frequency call of interest:
bandwidth_call_freq=[min tonal frequency; max tonal frequency]

window_call_durations=[min tonal duraitn; max tonal duration]
space_btw_diff_call= duration in sec between two energetic events to consider new call

thresh_sigma= best to set it at 2 or 3,
for instance thresh_sigma= 3 means  [−3 σ; 3 σ] with σ being the standard deviation of the distribution
here, for a random variable that follows a normal distribution, 99.74% of the values (corresponding
to ambient noise within our framework) belong to the [−3 σ; 3 σ] interval in comparison with the mean

appli_fit_type : can be set to 'no' to apply only the energetic detector,
otehrwise can be set to 'yes' or left empty
fit_type = 'poly3','exp1', 'exp2'... see "help fit"
thresh_Rsquare: minimum R² to consider fit ok

'threshold_detect' is estimated usin the mode, but can be change to median or mean if
requires, but mode seemed to be the best approach.


