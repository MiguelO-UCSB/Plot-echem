# Plot-echem

A Python graphical user interface designed to plot electrochemistry files for the Sepunaru lab.

## Installation

Clone the repo with GitHub.

```
  git clone https://github.com/MiguelO-UCSB/Plot-echem

```
    
## Usage/Examples

Open main.py in your favorite Python environment (Example: Spyder) or run with Anaconda Command Prompt.

Run the program in the environment.

![App Screenshot](https://github.com/MiguelO-UCSB/Plot-echem/blob/main/Other/GUI_example.png?raw=true)

Supported files:

   Biologic .txt
      
      - Format: Tab separated (time/s	Ewe/V	<I>/mA cycle number) or (time/s	Ewe/V	<I>/mA) or (time/s	Ewe/V) or (freq/Hz	Re(Z)/Ohm	-Im(Z)/Ohm	|Z|/Ohm	Phase(Z)/deg)
   
   SECCM .csv

      - Format: Comma separated (t/s,	E/V, I/A)
   
   Autolab EIS files .txt

      - Format: ; separated (Index;Frequency (Hz);Z' (Ω);-Z'' (Ω);Z (Ω);-Phase (°);Time (s)) or (Z' (Ω);-Z'' (Ω);Error Z' (%);Error Z'' (%);Z (Ω);-Phase (°);Frequency (Hz);Convergence;Number of iterations;χ²)
   
   HEKA .asc and .mat

      - ASC Format: Comma separated ("Index", "Time[s]", "Imon-1[A]", "Time[s]", "Vmon-1[V]")
