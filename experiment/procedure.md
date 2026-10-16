### Procedure

Follow these 8 experimental steps in the 3D Virtual Reality Simulation to analyze and mitigate 50 Hz Power Line Interference (PLI):

1. **Step 1: Configure the Bio-Signal Source**
   - In the controls dock, select your target physiological signal: **ECG** (electrocardiogram) or **EEG** (alpha rhythm).
   - For ECG, adjust the heart rate control (e.g. 72 BPM).
   - Note the baseline physiological amplitude ($1.0\text{ mV}$ peak-to-peak for ECG Lead II).

2. **Step 2: Generate & Observe the Clean Signal**
   - With power-line interference turned off ($0.0\text{ mV}$), observe Channel 1 (Cyan trace) on the Digital Storage Oscilloscope (DSO).
   - Verify that the cardiac P-wave, QRS complex, and T-wave are pristine with zero baseline hum.
   - Note the baseline Signal-to-Noise Ratio ($\text{SNR} > 35\text{ dB}$).

3. **Step 3: Introduce 50 Hz Power-Line Interference & Ground Loop**
   - Set the noise amplitude slider to $1.20\text{ mV}$.
   - Switch OFF Galvanic Isolation to simulate common-mode ground loop circulating currents.
   - Observe Channel 2 (Yellow trace) on the DSO. Notice how the ECG signal is completely corrupted and buried beneath the 50 Hz sinusoidal hum, driving the Input SNR to a negative value ($\text{SNR}_{\text{in}} < -4\text{ dB}$).

4. **Step 4: Inspect the Real-Time FFT Spectrum**
   - Navigate camera to the **FFT Spectrum Analyzer** (or click "Expanded DSO & FFT" in the top header).
   - Inspect the frequency distribution between $0\text{ Hz}$ and $150\text{ Hz}$.
   - Identify the tall fundamental power-line spike at exactly $50.0\text{ Hz}$ and the 3rd harmonic spike at $150.0\text{ Hz}$ towering over the physiological spectrum.

5. **Step 5: Enable & Configure the Digital IIR Notch Filter**
   - Toggle the **Active Filter State** to ON.
   - Select the 2nd-order **IIR Biquad** topology.
   - Set the notch center frequency $f_0 = 50.0\text{ Hz}$ to match the AC power grid frequency.
   - Adjust the Quality Factor $Q$ from $4.0$ to $16.0$. Observe the narrowing of the notch bandwidth ($\Delta f = f_0 / Q = 3.1\text{ Hz}$).
   - Re-engage **Optocoupler Galvanic Isolation** to break the ground loop and elevate the system CMRR to $> 120\text{ dB}$.

6. **Step 6: Compare Input and Filtered Signals on DSO**
   - Return to the Oscilloscope view and observe the 3 channels simultaneously:
     - **CH1 (Cyan):** Clean Reference Signal
     - **CH2 (Yellow):** Corrupted Noisy Input
     - **CH3 (Green):** Filtered Output Signal
   - Confirm that the notch filter eliminates the 50 Hz ripple while preserving the sharp R-wave amplitude and delicate P/T wave morphology.

7. **Step 7: Evaluate Quantitative Performance Metrics**
   - Check the live Performance Analytics HUD on the left:
     - **Input Noise RMS vs Output Error RMS** (e.g. $848.5\ \mu\text{V} \rightarrow 42.1\ \mu\text{V}$)
     - **Output SNR Improvement** ($\text{SNR}_{\text{out}} > +25\text{ dB}$)
     - **Estimated Noise Reduction** ($> 90\%$)
     - **50 Hz Notch Attenuation** ($> 45\text{ dB}$)
   - Click "Listen to 50Hz Hum" to verify audio silencing when the filter is active.

8. **Step 8: Record Observations & Export Results**
   - In Step 8 of the Guided Experiment wizard, review the comprehensive results summary table.
   - Enter your analytical notes regarding notch selectivity, $Q$-factor trade-offs, and patient electrical safety per IEC 60601-1.
   - Click **"Save CSV"** to download time-series numerical data or click **"Print / PDF Report"** to export a formal laboratory experiment report.