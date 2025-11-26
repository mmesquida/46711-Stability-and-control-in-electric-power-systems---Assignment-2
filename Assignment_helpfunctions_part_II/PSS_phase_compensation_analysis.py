import scipy.io as sio
import numpy as np
import control.statesp as stsp
import matplotlib.pyplot as plt
import stsp_functions as f_stsp
import os

""" Parameter definitions """
fmin = 0; fmax = 2
f = np.arange(fmin, fmax, .001)
omega = 2 * np.pi * f

# PSS design parameters
gen = 4  # Specify Generator where the PSS should be tuned
Ks = 20; Tw = 16; Tn1 = 0.0847; Td1 = 0.0000000001; Tn2 = 3.474; Td2 = 5.45

""" Load System Data """
current_path = os.getcwd()
parent_path = os.path.dirname(current_path)
file_path = os.path.join(parent_path, 'Assignment_data', 'system_q2.mat')

sys_data = sio.loadmat(file_path, squeeze_me=True)
A = sys_data['A']
B = sys_data['B']
C = sys_data['C']
D = sys_data['D']
strsps = sys_data['strsps']
statename = sys_data['StateName']

"""###################### Your code starts here ############################"""

""" Create PSS state space object """
pss = f_stsp.pss_stsp(Ks, Tw, Tn1, Td1, Tn2, Td2)

""" Find the location of relevant indices in the original system """
# Convert statename to string list for easier handling
statename_list = [str(name) for name in np.atleast_1d(statename)]

# Identify rotor states to exclude (angle and speed)
is_rotor = np.array([(name == 'angle') or (name == 'speed') for name in statename_list])
is_elec = ~is_rotor  # All non-rotor states

# Get input/output indices for specified generator
vref_idx = int(strsps['vref'].item()[gen-1] - 1)  # Vref input for generator
pe_idx = int(strsps['pe'].item()[gen-1] - 1)      # Pe output for generator

""" Manipulate the system by removing the dynamic states """
# Create reduced system matrices (excluding rotor states)
A_red = A[np.ix_(is_elec, is_elec)]
B_red = B[np.ix_(is_elec, [vref_idx])]
C_red = C[np.ix_([pe_idx], is_elec)]
D_red = D[np.ix_([pe_idx], [vref_idx])]

# Create reduced state space system
sys_red = stsp.StateSpace(A_red, B_red, C_red, D_red)

""" Determine the frequency response of the system and the PSS """
# System frequency response
mag_sys, phase_sys, _ = sys_red.frequency_response(omega, squeeze=False)
phase_sys_deg = np.degrees(np.squeeze(phase_sys))

# PSS frequency response  
mag_pss, phase_pss, _ = pss.frequency_response(omega, squeeze=False)
phase_pss_deg = np.degrees(np.squeeze(phase_pss))

# Wrap system phase to [-180, 180] range
phase_sys_deg = (phase_sys_deg + 180) % 360 - 180

# Ideal phase to be supplied = -system_phase (inverted as suggested by teacher)
ideal_phase = -phase_sys_deg

# Total phase after compensation
total_phase = phase_pss_deg + phase_sys_deg

""" Plot the frequency response curves """
plt.figure(figsize=(12, 8))

# Plot 1: Ideal phase compensation needed vs PSS provided phase
plt.subplot(2, 1, 1)
plt.plot(f, ideal_phase, 'b-', linewidth=2, label='Ideal Phase Compensation Needed')
plt.plot(f, phase_pss_deg, 'r--', linewidth=2, label='PSS Provided Phase')
plt.xlabel('Frequency [Hz]')
plt.ylabel('Phase [°]')
plt.title(f'PSS Phase Compensation Analysis - Generator {gen}')
plt.grid(True, alpha=0.3)
plt.legend()
plt.xlim(0, 2)
plt.ylim(bottom=0)
plt.tight_layout()
plt.show()

""" Analysis of compensation quality """
print(f"\nPSS Compensation Analysis for Generator {gen}")
print("="*50)
print(f"PSS Parameters: Ks={Ks}, Tw={Tw}, Tn1={Tn1}, Td1={Td1}, Tn2={Tn2}, Td2={Td2}")

key_frequencies = [0.3, 0.5, 0.8, 1.0, 1.2]
print("\nFrequency | Ideal Phase | PSS Phase | Error | Status")
print("-"*55)

for freq in key_frequencies:
    idx = np.argmin(np.abs(f - freq))
    ideal = ideal_phase[idx]
    pss_phase = phase_pss_deg[idx]
    error = total_phase[idx]
    
    if abs(error) <= 10:
        status = "GOOD"
    elif error > 10:
        status = "OVER-COMPENSATED"
    else:
        status = "UNDER-COMPENSATED"
    
    print(f"{freq:4.1f} Hz  | {ideal:8.1f}°    | {pss_phase:6.1f}°   | {error:5.1f}° | {status}")
