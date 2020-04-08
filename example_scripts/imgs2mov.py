# short script to generate a movie from plots in plot folder

import ffmpeg

(
    ffmpeg
    .input('plots/*.png', pattern_type='glob', framerate=3)
    .output('26kHz.mp4')
    .run()
)