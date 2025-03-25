# Differential-pulse-code-modulation
This project covers Differential Pulse Code Modulation (DPCM), a signal encoding method that predicts sample values based on previously encoded samples and transmits only the prediction error. The project involves:

    Implementing a DPCM encoder/decoder in MATLAB.

    Quantizing prediction errors using a uniform quantizer.

    Processing a 20,000-sample signal (source.mat) with different predictor orders and bit depths (N = 1,2,3).

    Evaluating performance through Mean Squared Error (MSE), signal reconstruction, and predictor coefficient analysis.

The goal is to analyze DPCM efficiency in reducing quantization error compared to standard PCM.
