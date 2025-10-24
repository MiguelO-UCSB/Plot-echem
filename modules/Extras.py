from tkinter import *
from tkinter.ttk import *
from tkinter import filedialog
from PIL import Image, ImageTk
import matplotlib
import os
from matplotlib import pyplot as plt, cm, ticker as tk
from matplotlib.backends.backend_tkagg import FigureCanvasTkAgg
from scipy import signal
from scipy.signal import find_peaks
import numpy as np

bool_map = {"true": True, 'True': True,
            "false": False, 'False': False}

Overlay_options = ['True', 'False']
Export_format_types = ['png', 'pdf', 'svg', 'jpg']

def OptionMenuStringVar(frame, options,row,col,sticky,idx=0, return_widget=False):
    var = StringVar()
    menu = OptionMenu(frame, var, options[idx],*options)
    menu.grid(row=row, column=col, sticky=sticky)
    if return_widget:
        return var, menu
    return var

class Make_Popup_GUI_Figure():
    
    def __init__(self, GUI):
        self.GUI = GUI
        self.fig2 = plt.Figure(figsize=(4,4), dpi=100, constrained_layout=True)
        # self.ax_signal = self.fig2.add_subplot(211)
        self.ax_fft = self.fig2.add_subplot(111)
        self.make_popup()
        
    def make_popup(self):
        self.popup = Toplevel()
        self.popup.title("Fourier Transform Data")
        self.leftframe  = Frame(self.popup)
        self.rightframe = Frame(self.popup)
        self.leftframe.grid(row=0, column=0)
        self.rightframe.grid(row=0, column=1)
        self.canvas = FigureCanvasTkAgg(self.fig2, master=self.rightframe)
        self.canvas.get_tk_widget().grid(row=0, column=0)
    
    def draw(self):
        pass
        
    def save(self):
        export_format = self.Export_format.get()
        path = filedialog.asksaveasfilename(defaultextension=export_format)
        if path == '':
            print('Error: plot not saved')
            return
        
        transparency = bool_map.get(self.Transparency.get().strip().lower(), False)
        dpi_val = int(self.dpi_field.get())
        
        self.fig2.savefig(path, transparent = transparency, dpi=dpi_val, format=export_format)
        print(f'Plot exported as {export_format}')

class FTdataPopup(Make_Popup_GUI_Figure):
    
    def __init__(self, GUI, echem_data):
        self.data = echem_data
        super().__init__(GUI=GUI)
        self.draw()
    
    def fill_leftframe(self):
        # Put relevant buttons in self.leftframe
        frame = self.leftframe
        
        if not hasattr(self, 'dpi_field'):
            # Will already have these attributes if it's being reinitialized
            self.prominence = StringVar(value='0.1')
            self.distance = StringVar(value='10')
            self.dpi_field = StringVar(value='300')
            self.xlabel = StringVar(value='Frequency (Hz)')
            self.ylabel = StringVar(value='Amplitude')
            self.x_axis_tickmultiple = StringVar(value='')
            self.x_axis_minor_tickmultiple = StringVar(value='')
            self.echem_xminval = StringVar(value='')
            self.echem_xmaxval = StringVar(value='')
            self.y_axis_tickmultiple = StringVar(value='')
            self.y_axis_minor_tickmultiple = StringVar(value='')
            self.echem_yminval = StringVar(value='')
            self.echem_ymaxval = StringVar(value='')
            self.n_yticks = StringVar(value='')
            self.div_const = StringVar(value='1')
        
        '''Make Frame for Voltamperometric Tab'''
        
        NormFrame = Frame(frame)
        NormFrame.grid(row=0, column=0, sticky=(N,S,E,W), pady=10)
        
        '''Norm'''
        Button(NormFrame, text='Redraw', command=self.redraw).grid(row=1, column=0, sticky=(W,E))
        Button(NormFrame, text='Save', command=self.save).grid(row=1, column=1, sticky=(W,E))
        Label(NormFrame, text='\nFind Peak Parameters').grid(row=2, column=0, sticky=(W,E), columnspan=2)
        Label(NormFrame, text='Prominence: ').grid(row=3, column=0, sticky=(E))
        Entry(NormFrame, textvariable=self.prominence, width=10).grid(
            row=3, column=1, columnspan=1, sticky=(W,E))
        Label(NormFrame, text='Distance: ').grid(row=4, column=0, sticky=(E))
        Entry(NormFrame, textvariable=self.distance, width=10).grid(
            row=4, column=1, columnspan=1, sticky=(W,E))
        
        
        # Make Nested Notebooks
        inner_tabs = Notebook(frame)
        inner_tabs.grid(row=1, column=0, sticky=(N,S,E,W), pady=10)
        frame.grid_rowconfigure(0, weight=1)
        frame.grid_columnconfigure(0, weight=1)
        AxisFrame = Frame(inner_tabs)
        ExportFrame = Frame (inner_tabs)
        
        inner_tabs.add(AxisFrame, text='  Axis  ')
        inner_tabs.add(ExportFrame, text='  Export  ')
        
        '''Axis'''
        Label(AxisFrame, text='X-axis label: ').grid(row=0, column=0, columnspan=1, sticky=(E))
        Entry(AxisFrame, textvariable=self.xlabel, width=10).grid(
            row=0, column=1, columnspan=1, sticky=(W,E))
        Label(AxisFrame, text='Y-axis label: ').grid(row=1, column=0, columnspan=1, sticky=(E))
        Entry(AxisFrame, textvariable=self.ylabel, width=10).grid(
            row=1, column=1, columnspan=1, sticky=(W,E))
        Label(AxisFrame, text='X-axis tick multiple: ').grid(row=2, column=0, sticky=(E))
        Entry(AxisFrame, textvariable=self.x_axis_tickmultiple, width=10).grid(
            row=2, column=1, columnspan=1, sticky=(W,E))
        Label(AxisFrame, text='Minor tick multiple: ').grid(row=3, column=0, sticky=(E))
        Entry(AxisFrame, textvariable=self.x_axis_minor_tickmultiple, width=10).grid(
            row=3, column=1, columnspan=1, sticky=(W,E))
        Label(AxisFrame, text='Min: ').grid(row=4, column=0, sticky=(E))
        Entry(AxisFrame, textvariable=self.echem_xminval, width=10).grid(
            row=4, column=1, columnspan=1, sticky=(W,E))
        Label(AxisFrame, text='Max: ').grid(row=4, column=3, sticky=(E))
        Entry(AxisFrame, textvariable=self.echem_xmaxval, width=10).grid(
            row=4, column=4, columnspan=1, sticky=(W,E))
        
        Label(AxisFrame, text='').grid(row=5, column=0, sticky=(E))
        
        Label(AxisFrame, text='Y-axis tick multiple: ').grid(row=9, column=0, sticky=(E))
        Entry(AxisFrame, textvariable=self.y_axis_tickmultiple, width=10).grid(
            row=9, column=1, columnspan=1, sticky=(W,E))
        Label(AxisFrame, text='Minor tick multiple: ').grid(row=10, column=0, sticky=(E))
        Entry(AxisFrame, textvariable=self.y_axis_minor_tickmultiple, width=10).grid(
            row=10, column=1, columnspan=1, sticky=(W,E))
        Label(AxisFrame, text='Min: ').grid(row=11, column=0, sticky=(E))
        Entry(AxisFrame, textvariable=self.echem_yminval, width=10).grid(
            row=11, column=1, columnspan=1, sticky=(W,E))
        Label(AxisFrame, text='Max: ').grid(row=11, column=3, sticky=(E))
        Entry(AxisFrame, textvariable=self.echem_ymaxval, width=10).grid(
            row=11, column=4, columnspan=1, sticky=(W,E))
        Label(AxisFrame, text='').grid(row=12, column=0, sticky=(E))
        Label(AxisFrame, text='Y-axis divide constant: ').grid(row=13, column=0, 
                                                    columnspan=1, sticky=(E))
        Entry(AxisFrame, textvariable=self.div_const, width=10).grid(
            row=13, column=1, columnspan=2, sticky=(W,E))
        Label(AxisFrame, text='').grid(row=14, column=0, sticky=(E))
        
        '''Export'''
        #Trans - True/False
        #dpi #
        #format - 'png', 'pdf', 'svg', 'jpg'
        Label(ExportFrame, text='Transparent: ').grid(row=0, column=0, sticky=(E))
        self.Transparency = OptionMenuStringVar(ExportFrame, Overlay_options, 0, 1, (W,E), idx=1,)
        
        Label(ExportFrame, text='dpi: ').grid(row=1, column=0, sticky=(E))
        Entry(ExportFrame, textvariable=self.dpi_field, width=4).grid(row=1, column=1,
                                                                sticky=(W,E))
        
        Label(ExportFrame, text='Format: ').grid(row=3, column=0, sticky=(E))
        self.Export_format = OptionMenuStringVar(ExportFrame, Export_format_types, 3, 1, (W,E), idx=0)
        
        Label(frame, text='').grid(row=4, column=0, sticky=(E))

        pass
    
    def redraw(self):
        print('Fourier Transform Plot Redrawn')
        self.ax_fft.clear() 
        self.draw()
    
    
    def draw(self):
        # ---- Step 1: Combine signals into one (for FT) ----
        # Flatten all x, y arrays into one long signal
        x_combined = np.concatenate(self.data["x"])
        y_combined = np.concatenate(self.data["y"])

        # Ensure sorted by x (important for FFT spacing)
        sort_idx = np.argsort(x_combined)
        x_sorted = x_combined[sort_idx]
        y_sorted = y_combined[sort_idx]

        # ---- Step 2: Do Fourier Transform ----
        # Assume uniform spacing
        dx = np.mean(np.diff(x_sorted))
        N = len(y_sorted)

        freqs = np.fft.fftfreq(N, d=dx)
        fft_vals = np.fft.fft(y_sorted)

        # Only keep positive frequencies
        mask = freqs >= 0
        freqs = freqs[mask]
        fft_vals = np.abs(fft_vals[mask])

        # ---- Step 3: Plot result in popup ----
        # FT Magnitude
        self.ax_fft.plot(freqs, fft_vals/float(self.div_const.get()), lw=2, color='k')
        
        self.find_peaks(freqs, fft_vals/float(self.div_const.get()))
        
        self.set_xlabel()
        self.set_ylabel()
        self.update_plotlim()
        
        self.canvas.draw()
    
    def set_xlabel(self):
        xlabel = self.xlabel.get()
        self.ax_fft.set_xlabel(f'{xlabel}')
    
    def set_ylabel(self):
        ylabel = self.ylabel.get()
        self.ax_fft.set_ylabel(f'{ylabel}')
    
    def update_plotlim(self):
        '''
        Set Axis limits and tick multiples for IV plots
        '''
        # Grab values from GUI ('' if empty)
        x_multiple = self.x_axis_tickmultiple.get()
        xmin = self.echem_xminval.get()
        xmax = self.echem_xmaxval.get()
        
        y_multiple = self.y_axis_tickmultiple.get()
        ymin = self.echem_yminval.get()
        ymax = self.echem_ymaxval.get()
        # lim = (xmin, xmax), (ymin, ymax)
        
        x_minor_multiple = self.x_axis_minor_tickmultiple.get()
        y_minor_multiple = self.y_axis_minor_tickmultiple.get()
        
        xlim = (xmin, xmax)
        ylim = (ymin, ymax)
        
        # Update x-axis limits if both are provided
        if all(x != '' for x in xlim):
            try:
                self.ax_fft.set_xlim(float(xlim[0]), float(xlim[1]))
            except Exception as e:
                print(f'Skiping setting x axis limit because of error: {e}')
                pass  # Could log or print warning
        
        # Update y-axis limits if both are provided
        if all(y != '' for y in ylim):
            try:
                self.ax_fft.set_ylim(float(ylim[0]), float(ylim[1]))
            except Exception as e:
                print(f'Skiping setting y axis limit because of error: {e}')
                pass
        
        # Update x-axis tick multiple if provided
        if x_multiple != '':
            try:
                step = float(x_multiple)
                xmin, xmax = self.ax_fft.get_xlim()
                n_ticks = abs((xmax - xmin) / step)
                if n_ticks > 100:
                    print(f"Error ⚠️ Too many major x-ticks ({int(n_ticks)}). Using AutoLocator instead...")
                    self.ax_fft.xaxis.set_major_locator(tk.AutoLocator())
                else:
                    self.ax_fft.xaxis.set_major_locator(tk.MultipleLocator(float(x_multiple)))
            except Exception as e:
                print(f'Skiping setting x axis major tick multiple because of error: {e}')
                pass
        
        # Update y-axis tick multiple if provided
        if y_multiple != '':
            try:
                step = float(y_multiple)
                ymin, ymax = self.ax_fft.get_ylim()
                n_ticks = abs((ymax - ymin) / step)
                if n_ticks > 100:
                    print(f"Error ⚠️ Too many major y-ticks ({int(n_ticks)}). Using AutoLocator instead...")
                    self.ax_fft.yaxis.set_major_locator(tk.AutoLocator())
                else:
                    self.ax_fft.yaxis.set_major_locator(tk.MultipleLocator(float(y_multiple)))
            except Exception as e:
                print(f'Skiping setting y axis major tick multiple because of error: {e}')
                pass
        
        # Update x-axis minor tick multiple if provided
        if x_minor_multiple != '':
            try:
                step = float(x_minor_multiple)
                xmin, xmax = self.ax_fft.get_xlim()
                n_ticks = abs((xmax - xmin) / step)
                if n_ticks > 100:
                    print(f"Error ⚠️ Too many minor x-ticks ({int(n_ticks)}). Using AutoLocator instead...")
                    self.ax_fft.xaxis.set_minor_locator(tk.AutoLocator())
                else:
                    self.ax_fft.xaxis.set_minor_locator(tk.MultipleLocator(float(x_minor_multiple)))
            except Exception as e:
                print(f'Skiping setting x axis minor tick multiple because of error: {e}')
                pass
        
        # Update y-axis minor tick multiple if provided
        if y_minor_multiple != '':
            try:
                step = float(y_minor_multiple)
                ymin, ymax = self.ax_fft.get_ylim()
                n_ticks = abs((ymax - ymin) / step)
                if n_ticks > 100:
                    print(f"Error ⚠️ Too many minor y-ticks ({int(n_ticks)}). Using AutoLocator instead...")
                    self.ax_fft.yaxis.set_minor_locator(tk.AutoLocator())
                else:
                    self.ax_fft.yaxis.set_minor_locator(tk.MultipleLocator(float(y_minor_multiple)))
            except Exception as e:
                print(f'Skiping setting y axis minor tick multiple because of error: {e}')
                pass
    
    def find_peaks(self, freqs_fil, ft_V_fil):
        '''Find peak locations'''
        try:
            prominence_ = float(self.prominence.get())
            distance_ = float(self.distance.get())
            peaks1, props1 = find_peaks(abs(ft_V_fil), prominence=prominence_, distance=distance_, )
        except Exception as e:
            print(f'Error: Could not find peak locations because of error {e}')
            return
        
        if len(peaks1) == 0:
            print('No peaks found. Change prominence')
            return
        
        '''Filter duplicates'''
        deduped = []
        min_step = 10
        # x_threshold = 1
        x_threshold = False
        current_group = [peaks1[0]]
        for i in range(1, len(peaks1)):
            if freqs_fil[peaks1[i]] - freqs_fil[current_group[-1]] <= min_step:
                current_group.append(peaks1[i])
            else:
                # Choose the peak with the highest y-value
                max_peak = max(current_group, key=lambda i: abs(ft_V_fil)[i])
                if x_threshold is False or freqs_fil[max_peak] >= x_threshold:
                    deduped.append(max_peak)
                current_group = [peaks1[i]]
        
        if current_group: # Don't forget the last group
            max_peak = max(current_group, key=lambda i: abs(ft_V_fil)[i])
            if x_threshold is False or freqs_fil[max_peak] >= x_threshold:
                deduped.append(max_peak)
        
        self.ax_fft.plot(freqs_fil[deduped], abs(ft_V_fil)[deduped],
                color = 'r',
                marker='o', markersize=3,
                linewidth = 0)
        
        for i in range(len(freqs_fil[deduped][freqs_fil[deduped]<60])):
            self.ax_fft.text(freqs_fil[deduped][i]+1, abs(ft_V_fil)[deduped][i],
                             f'{round(freqs_fil[deduped][i], 2)} Hz',
                             ha='left', size=12, rotation='horizontal')
    
class Make_Popup_GUI():
    
    def __init__(self, GUI):
        self.GUI = GUI
        self.make_popup()
        
        # --- Data: Potentials vs. NHE at 25 °C in Volts ---
        # These are the standard values used for conversion.
        self.REF_POTENTIALS = {
            "NHE (Normal Hydrogen Electrode)": 0.0000,
            "SHE (Standard Hydrogen Electrode)": 0.0000,
            "SCE (Saturated Calomel Electrode)": 0.2412,
            "Calomel (0.1 M KCl)": 0.3337,
            "Calomel (1.0 M KCl)": 0.2801,
            "Ag/AgCl (Saturated KCl)": 0.197,
            "Ag/AgCl (0.1 M KCl)": 0.2881,
            "Ag/AgCl (1.0 M KCl)": 0.235,
            "Hg/Hg₂SO₄ (Saturated K₂SO₄)": 0.640,
            "Hg/HgO (0.1 M NaOH)": 0.149,
            "Hg/HgO (1.0 M NaOH)": 0.108,
            "Fc/Fc⁺ (Ferrocene/Ferrocenium in MeCN)": 0.640
        }
        
        self.REF_POTENTIALS_RHE = {
            "SHE (Standard Hydrogen Electrode)": 0,
            "RHE (Reversible Hydrogen Electrode)": 0}

        # --- Tkinter String Variables ---
        # These variables are linked to widgets to get/set their values easily.
        if not hasattr(self, 'input_potential_var'):
            # Will already have these attributes if it's being reinitialized
            self.input_potential_var = StringVar(value='0')
            self.input_potential_var_RHE = StringVar(value='0')
            self.temperature = StringVar(value='298')
            self.pH = StringVar(value='1')
            
        self.input_ref_var = StringVar()
        self.input_ref_var_RHE = StringVar()
        
        # Create a dictionary of StringVars for each output label
        self.output_vars = {name: StringVar(value="--.--") for name in self.REF_POTENTIALS}
        self.output_vars_RHE = {name: StringVar(value="--.--") for name in self.REF_POTENTIALS_RHE}

        # Set a default value for the dropdown menu
        self.input_ref_var.set(list(self.REF_POTENTIALS.keys())[2]) # Default to SCE
        self.input_ref_var_RHE.set(list(self.REF_POTENTIALS_RHE.keys())[0]) # Default to SCE
        
        self.fill_leftframe()
        self.fill_rightframe()
        
    def make_popup(self):
        self.popup = Toplevel()
        self.popup.title("Reference Electrode Converter")
        self.leftframe  = Frame(self.popup)
        self.rightframe = Frame(self.popup)
        self.leftframe.grid(row=0, column=0)
        self.rightframe.grid(row=0, column=1)
        
        
    def fill_leftframe(self):
        pass
    
    def fill_rightframe(self):
        pass


class ReferenceElecConvPopup(Make_Popup_GUI):
    
    def __init__(self, GUI):
        super().__init__(GUI=GUI)
    
    def fill_leftframe(self):
        # Put relevant buttons in self.leftframe
        frame = self.leftframe
        
        # --- Input Section ---
        input_frame = LabelFrame(frame, text="Input", padding="10")
        input_frame.grid(row=0, column=0, sticky="ew", padx=5, pady=5)
        input_frame.grid_columnconfigure(1, weight=1)

        Label(input_frame, text="Potential (V):").grid(row=0, column=0, padx=5, pady=5, sticky="w")
        potential_entry = Entry(input_frame, textvariable=self.input_potential_var)
        potential_entry.grid(row=0, column=1, padx=5, pady=5, sticky="ew")

        Label(input_frame, text="Reference:").grid(row=1, column=0, padx=5, pady=5, sticky="w")
        ref_dropdown = OptionMenu(input_frame, self.input_ref_var, self.input_ref_var.get(), *self.REF_POTENTIALS.keys())
        ref_dropdown.grid(row=1, column=1, padx=5, pady=5, sticky="ew")

        # --- Convert Button ---
        convert_button = Button(frame, text="Convert", command=self.convert_potentials)
        convert_button.grid(row=1, column=0, padx=15, pady=10, sticky="ew")

        # --- Output Section ---
        output_frame = LabelFrame(frame, text="Converted Potentials", padding="10")
        output_frame.grid(row=2, column=0, sticky="nsew", padx=5, pady=5)
        output_frame.grid_columnconfigure(1, weight=1)
        
        # Dynamically create labels for each reference electrode
        for i, name in enumerate(self.REF_POTENTIALS.keys()):
            Label(output_frame, text=f"{name}:").grid(row=i, column=0, padx=5, pady=2, sticky="w")
            Label(output_frame, textvariable=self.output_vars[name], font=("Calibri", 14)).grid(row=i, column=1, padx=5, pady=2, sticky="w")
        return
        
    def convert_potentials(self):
        """The core logic for the conversion."""
        try:
            # 1. Get the user's input potential as a number
            input_potential = float(self.input_potential_var.get())
            
            # 2. Get the name of the user's reference electrode
            input_ref_name = self.input_ref_var.get()
            
            # 3. Get the potential of the user's reference vs. NHE
            ref_potential_vs_nhe = self.REF_POTENTIALS[input_ref_name]

            # 4. Convert the input potential to the NHE scale
            # E_vs_NHE = E_vs_Ref + E_Ref_vs_NHE
            potential_vs_nhe = input_potential + ref_potential_vs_nhe

            # 5. Convert from NHE to all other reference scales and update the labels
            # E_vs_NewRef = E_vs_NHE - E_NewRef_vs_NHE
            for name, value_vs_nhe in self.REF_POTENTIALS.items():
                converted_value = potential_vs_nhe - value_vs_nhe
                self.output_vars[name].set(f"{converted_value:.4f} V") # Format to 4 decimal places
        
        except ValueError:
            # Handle cases where the input is not a valid number
            for name in self.REF_POTENTIALS:
                self.output_vars[name].set("Invalid Input")
        except Exception as e:
            # Handle any other unexpected errors
            print(f"An error occurred: {e}")
            for name in self.REF_POTENTIALS:
                self.output_vars[name].set("Error")
    
    def fill_rightframe(self):
        # Put relevant buttons in self.leftframe
        frame = self.rightframe
        
        # --- Input Section ---
        input_frame = LabelFrame(frame, text="Input", padding="10")
        input_frame.grid(row=0, column=0, sticky="ew", padx=5, pady=5)
        input_frame.grid_columnconfigure(1, weight=1)

        Label(input_frame, text="Potential (V):").grid(row=0, column=0, padx=5, pady=5, sticky="w")
        potential_entry = Entry(input_frame, textvariable=self.input_potential_var_RHE)
        potential_entry.grid(row=0, column=1, padx=5, pady=5, sticky="ew")

        Label(input_frame, text="Reference:").grid(row=1, column=0, padx=5, pady=5, sticky="w")
        ref_dropdown = OptionMenu(input_frame, self.input_ref_var_RHE, self.input_ref_var_RHE.get(), *self.REF_POTENTIALS_RHE.keys())
        ref_dropdown.grid(row=1, column=1, padx=5, pady=5, sticky="ew")
        
        Label(input_frame, text="Temperature (K):").grid(row=2, column=0, padx=5, pady=5, sticky="w")
        potential_entry = Entry(input_frame, textvariable=self.temperature)
        potential_entry.grid(row=2, column=1, padx=5, pady=5, sticky="ew")
        
        Label(input_frame, text="pH:").grid(row=3, column=0, padx=5, pady=5, sticky="w")
        potential_entry = Entry(input_frame, textvariable=self.pH)
        potential_entry.grid(row=3, column=1, padx=5, pady=5, sticky="ew")

        # --- Convert Button ---
        convert_button = Button(frame, text="Convert", command=self.convert_RHE_potentials)
        convert_button.grid(row=1, column=0, padx=15, pady=10, sticky="ew")

        # --- Output Section ---
        output_frame = LabelFrame(frame, text="Converted Potential", padding="10")
        output_frame.grid(row=2, column=0, sticky="nsew", padx=5, pady=5)
        output_frame.grid_columnconfigure(1, weight=1)
        
        for i, name in enumerate(self.REF_POTENTIALS_RHE.keys()):
            Label(output_frame, text=f"{name}:").grid(row=i, column=0, padx=5, pady=2, sticky="w")
            Label(output_frame, textvariable=self.output_vars_RHE[name], font=("Calibri", 14)).grid(row=i, column=1, padx=5, pady=2, sticky="w")
        return
    
    def convert_RHE_potentials(self):
        """The core logic for the conversion."""
        try:
            # 1. Get the user's input potential as a number
            input_potential = float(self.input_potential_var_RHE.get())
            
            # 2. Get the name of the user's reference electrode and other params
            input_ref_name = self.input_ref_var_RHE.get()
            pH_ = float(self.pH.get())
            temp = float(self.temperature.get())
            
            if input_ref_name == "SHE (Standard Hydrogen Electrode)":
                converted_value = input_potential + np.log(10)*8.31446261815324*temp*pH_/96485.33128
                self.output_vars_RHE["SHE (Standard Hydrogen Electrode)"].set(f"{input_potential:.4f} V") # Format to 4 decimal places
                self.output_vars_RHE["RHE (Reversible Hydrogen Electrode)"].set(f"{converted_value:.4f} V") # Format to 4 decimal places
                
            if input_ref_name == "RHE (Reversible Hydrogen Electrode)":
                converted_value = input_potential - np.log(10)*8.31446261815324*temp*pH_/96485.33128
                self.output_vars_RHE["RHE (Reversible Hydrogen Electrode)"].set(f"{input_potential:.4f} V") # Format to 4 decimal places
                self.output_vars_RHE["SHE (Standard Hydrogen Electrode)"].set(f"{converted_value:.4f} V") # Format to 4 decimal places
                
        except ValueError:
            # Handle cases where the input is not a valid number
            for name in self.REF_POTENTIALS_RHE:
                self.output_vars_RHE[name].set("Invalid Input")
        except Exception as e:
            # Handle any other unexpected errors
            print(f"An error occurred: {e}")
            for name in self.REF_POTENTIALS_RHE:
                self.output_vars_RHE[name].set("Error")

class Make_Popup_GUI_macro():
    
    def __init__(self, GUI):
        self.GUI = GUI
        self.make_popup()
        
        # --- Tkinter String Variables ---
        # These variables are linked to widgets to get/set their values easily.
        if not hasattr(self, 'input_potential_var'):
            # Will already have these attributes if it's being reinitialized
            self.input_potential_var = StringVar(value='0')
        
        self.fill_leftframe()
        
    def make_popup(self):
        self.popup = Toplevel()
        self.popup.title("Macrodisk Randles–Ševčík Calculator")
        self.leftframe  = Frame(self.popup)
        self.rightframe = Frame(self.popup)
        self.leftframe.grid(row=0, column=0)
        self.rightframe.grid(row=0, column=1)
    
class MacroRSPopup(Make_Popup_GUI_macro):
    
    def __init__(self, GUI):
        super().__init__(GUI=GUI)
    
    def fill_leftframe(self):
        # Put relevant buttons in self.leftframe
        frame = self.leftframe
        
        # --- Reversible Image ---
        script_dir = os.path.dirname(os.path.abspath(__file__))
        macro_path = os.path.join(script_dir, "Images", "Rev_macro.png")
        macro_path_irr = os.path.join(script_dir, "Images", "Irr_macro.png")

        image = Image.open(macro_path)
        w, h = image.size
        image = image.resize((int(w * 0.32), int(h * 0.32)), Image.LANCZOS)  # 32% size
        image.load()  # force decode
        image = image.convert("RGBA")
        self.photo_image_rev = ImageTk.PhotoImage(image)
        
        image = Image.open(macro_path_irr)
        w, h = image.size
        image = image.resize((int(w * 0.32), int(h * 0.32)), Image.LANCZOS)  # 32% size
        image.load()  # force decode
        image = image.convert("RGBA")
        self.photo_image_irr = ImageTk.PhotoImage(image)
        
        Label(frame, text="Reversible", font=(8)).grid(row=1, column=0, padx=5, pady=5, sticky="w")
        Label(frame, image=self.photo_image_rev).grid(row=2, column=0, padx=5, pady=5, sticky="w")
        Label(frame, text="Irreversible", font=(8)).grid(row=3, column=0, padx=5, pady=5, sticky="w")
        Label(frame, image=self.photo_image_irr).grid(row=4, column=0, padx=5, pady=5, sticky="w")
        
        # --- Input Section ---
        input_frame = LabelFrame(frame, text="Input", padding="10")
        input_frame.grid(row=5, column=0, sticky="ew", padx=5, pady=5)
        input_frame.grid_columnconfigure(1, weight=1)
        
        return
    
class Popup_Generator():
    '''
    Class to run reference electrode converter
    '''
    
    def __init__(self):
        self.ReferenceElecConvPopup = None
        self.FTdataPopup = None
        self.MacroRSPopup = None
    
    def get(self, GUI, popup_type, data):
        '''
        Either reinitialize to remake window with old settings, or create a new object
        '''
        if popup_type == 'Ref_converter':
            if self.ReferenceElecConvPopup:
                # print('Test Reload')
                return self.ReferenceElecConvPopup.__init__(GUI)
            # print('Test First Open')
            self.ReferenceElecConvPopup = ReferenceElecConvPopup(GUI)
            return
        
        if popup_type == 'FT_data':
            if self.FTdataPopup:
                # print('Test Reload')
                return self.FTdataPopup.__init__(GUI, data)
            # print('Test First Open')
            self.FTdataPopup = FTdataPopup(GUI, data)
            return
        
        if popup_type == 'Macro_RS':
            if self.MacroRSPopup:
                # print('Test Reload')
                return self.MacroRSPopup.__init__(GUI)
            # print('Test First Open')
            self.MacroRSPopup = MacroRSPopup(GUI)
            return