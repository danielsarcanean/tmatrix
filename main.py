import numpy as np
import matplotlib.pyplot as plt
import csv

###----Constants----###

plt.rcParams["figure.dpi"] = 300
fig, ax = plt.subplots(figsize=(10,5))
me =  1 # SI: 9.11E-31 kg
hbar =  1 / (2 * np.pi) # SI: 6.626E-34 / (2 * np.pi) J · s
energy = 20 / 27.21 # --> convertion to hartrees from eV 
A = 0.1 # incident amplitude

###-----------------###

def tMatrixPositiveEnergy(xa, xb, amplitude):
    length = positions[xb] - positions[xa]
    A = amplitude
    if energy < potentials[xa]:
        k = np.sqrt(2 * me * (potentials[xa] - energy) / hbar**2)
        matrix = np.array([
            [0, 0],
            [0, np.exp(-k*length)]
        ])
        x = np.linspace(positions[xa], positions[xb+1], 200)
        plt.plot(x, energy + A*np.exp(-k*(x)), color='red')
    
    elif energy > potentials[xa] and energy < potentials[xb]:
        k1 = np.sqrt(2 * me * (energy - potentials[xa]) / hbar**2)
        k2 = np.sqrt(2 * me * (potentials[xb] - energy) / hbar**2)
        matrix = np.array([
            [0, (k1+1j*k2) / 2 * k1],
            [0, (k1-1j*k2) / 2 * k1]
        ])
        x = np.linspace(positions[xa], positions[xb+1], 200)
        plt.plot(x, energy + A*np.exp(-k2*(x - positions[xa])), color='red')
    
    elif potentials[xa] == potentials[xb]:
        k = np.sqrt(2 * me * (energy - potentials[xa]) / hbar**2)
        matrix = np.array([
            [np.exp(1j*k*length), 0],
            [0, np.exp(-1j*k*length)]
        ])
        x = np.linspace(positions[xa], positions[xb], 200)
        #plt.plot(x, energy +(energy - potentials[xb])*np.exp(-1j*k*(x)), color='green')
        plt.plot(x, energy + A*np.exp(1j*k*(x- positions[xa])), color='red')
    
    else:
        k1 = np.sqrt(2 * me * (energy - potentials[xa]) / hbar**2)
        k2 = np.sqrt(2 * me * (energy - potentials[xb]) / hbar**2)
        matrix = (0.5)*np.array([
            [1 + k2/k1, 1 - k2/k1],
            [1 - k2/k1, 1 + k2/k1]
        ])
    
    return matrix

def graphicSolution(potentials, positions):
    for i in range(len(potentials)-1):
        if potentials[i] == potentials[i+1]:
            ax.hlines(y=potentials[i], xmin=positions[i], xmax=positions[i+1])
        else: 
            ax.vlines(x=positions[i], ymin=potentials[i], ymax=potentials[i+1])
        ax.hlines(y=energy, xmin=positions[i], xmax=positions[i+1], linestyle='--', color='green')
    plt.ylabel("Potential (hartrees)")
    plt.xlabel("Position (bohr radius)")
    plt.legend(frameon=False, fontsize=9)
    plt.show()


def transferParameters(xa, xb, amplitude):
    initial_amplitude = amplitude
    finalMatrix = np.array([[1, 0], [0, 1]])
    for i in range(xb - xa):
        nextMatrix = tMatrixPositiveEnergy(i, i+1, amplitude)
        finalMatrix = np.matmul(nextMatrix, finalMatrix)
        amplitude = (np.linalg.det(finalMatrix) / finalMatrix[1,1]) * amplitude
        if energy < potentials[i]:
            break
        elif energy > potentials[i] and energy < potentials[i+1]:
            break
        else:
            continue
    transmitance = np.linalg.det(finalMatrix) / finalMatrix[1,1]
    reflectance = finalMatrix[1,0] / finalMatrix[1,1]
    print(f''' 
          Welcome to the t-Matrix calculator!
          
          ---------------------------------------------
          -            Your chosen output             -
          -         should have been plotted          -
          ---------------------------------------------
          
          Your initial amplitude is {initial_amplitude:.4f},
          and your transmitted and reflected amplitudes are 
          {amplitude:.4f} and {amplitude*reflectance:.4f} respectively.
          
          In terms of probabilities, your transmitance is {np.abs(transmitance)**2:.4f}
          and your reflectance is {np.abs(reflectance)**2:.4f}. Relating this to your
          amplitude's probability you have {np.abs(initial_amplitude)**2 * np.abs(transmitance)**2:.4f} for the
          Transmitance and {np.abs(initial_amplitude)**2 * np.abs(reflectance)**2:.4f} for the Reflectance.
          
          ''')
    return finalMatrix

### Variables of the program ### 

potentials = []
positions = []

with open('potentials.csv', 'r') as csv_file:
    csv_reader = csv.DictReader(csv_file, delimiter=';', skipinitialspace=True)
        
    for row in csv_reader:
        potentials.append(float(row['Potentials (eV)'])/27.21) # We ensure the eV values are passed to hartrees 
        positions.append(float(row['Positions (nm)'])*18.897) # in bohr radius

### Execution ### 

transferParameters(0,11,A)
graphicSolution(potentials, positions)
