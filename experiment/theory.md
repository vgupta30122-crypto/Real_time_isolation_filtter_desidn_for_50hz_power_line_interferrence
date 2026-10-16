### Theory

#### 1. Introduction & Origin of 50 Hz Power Line Interference (PLI)
In biopotential acquisition systems such as electrocardiography (ECG), electroencephalography (EEG), and electromyography (EMG), physiological signal amplitudes are exceptionally small:
- **ECG signals:** typically $0.5\text{ mV}$ to $2.0\text{ mV}$ peak-to-peak.
- **EEG signals:** typically $10\ \mu\text{V}$ to $100\ \mu\text{V}$.

In residential and clinical environments, 230 V AC / 50 Hz power distribution lines radiate time-varying electromagnetic fields. The human body acts as an antenna with an effective capacitance to the power line ($C_{\text{mains}} \approx 3\text{ pF} - 10\text{ pF}$) and capacitance to earth ground ($C_{\text{body}} \approx 100\text{ pF} - 250\text{ pF}$). This capacitive coupling induces displacement currents:

$$i_{\text{displacement}} = C_{\text{mains}} \cdot \frac{d V_{\text{mains}}}{dt} = C_{\text{mains}} \cdot 2\pi f_{\text{mains}} V_{\text{peak}}$$

These currents flow through electrode-skin impedances and ground loops, creating massive common-mode and differential 50 Hz voltages that can completely obscure low-amplitude diagnostic waveforms.

---

#### 2. Galvanic Isolation Barrier & Common Mode Rejection Ratio (CMRR)
A critical defense against mains interference and electric shock hazard is **Galvanic Isolation**.
- **Ground Loops:** When a patient recording circuit shares an earth ground with mains-powered recording instruments, potential differences between different ground points generate circulating loop currents, converting common-mode interference into large differential artifacts.
- **Galvanic Isolation Barrier:** An optocoupler (optical barrier) or capacitive isolation amplifier separates the patient-connected front end (floating ground) from the mains-grounded digital and acquisition system.
- **IEC 60601-1 Standard:** Provides $> 4\text{ kV}$ isolation dielectric withstand and reduces leakage currents below $10\ \mu\text{A}$.
- **CMRR Enhancement:**
  $$\text{CMRR (dB)} = 20 \log_{10}\left( \frac{A_d}{A_{cm}} \right) > 120\text{ dB}$$
  This effectively attenuates common-mode interference by over $1,000,000 : 1$.

---

#### 3. Real-Time 50 Hz Notch Filter Design

##### A. Digital 2nd-Order IIR Biquad Notch Filter
The discrete-time transfer function of a digital notch filter centered at $\omega_0 = 2\pi f_0 / f_s$ is given by:

$$H(z) = b_0 \frac{1 - 2\cos(\omega_0)z^{-1} + z^{-2}}{1 - 2r\cos(\omega_0)z^{-1} + r^2 z^{-2}}$$

Where:
- $f_0 = 50.0\text{ Hz}$ (Notch center frequency)
- $f_s = 1000\text{ Hz}$ (Biomedical sampling frequency)
- $\omega_0 = \frac{2\pi f_0}{f_s}$
- The 3-dB bandwidth $\Delta f = f_0 / Q$
- The pole radius $r \approx 1 - \frac{\pi \Delta f}{f_s}$

In Direct Form II Transposed implementation, the difference equations are:
$$y[n] = b_0 x[n] + w_1[n-1]$$
$$w_1[n] = b_1 x[n] - a_1 y[n] + w_2[n-1]$$
$$w_2[n] = b_2 x[n] - a_2 y[n]$$

##### B. Active Twin-T Analog Filter
The analog Twin-T filter consists of two parallel T-networks:
- A low-pass T-network: two resistors $R$ and a capacitor $2C$.
- A high-pass T-network: two capacitors $C$ and a resistor $R/2$.
The resonant notch frequency is:
$$f_0 = \frac{1}{2\pi R C} = 50.0\text{ Hz}$$
Bootstrapping with an operational amplifier (e.g., TL072 / OP07) increases the effective $Q$ from $0.25$ up to $> 15$, narrowing the notch so vital ECG frequency content remains uncorrupted.

##### C. Adaptive LMS Noise Cancellation
Using synthetic orthogonal reference vectors $x_1[n] = \cos(2\pi f_0 t)$ and $x_2[n] = \sin(2\pi f_0 t)$, the weights are updated iteratively:
$$e[n] = x_{\text{raw}}[n] - (w_1[n] x_1[n] + w_2[n] x_2[n])$$
$$w_k[n+1] = w_k[n] + 2\mu e[n] x_k[n]$$
This adaptively removes the 50 Hz sinusoid without introducing non-linear phase distortions into the QRS complex.