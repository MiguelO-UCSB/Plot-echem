from tkinter import *
from tkinter.ttk import *
from tkinter import filedialog
import matplotlib
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
        self.fill_leftframe()
        
    def make_popup(self):
        self.popup = Toplevel()
        self.popup.title("Fourier Transform Data")
        self.leftframe  = Frame(self.popup)
        self.rightframe = Frame(self.popup)
        self.leftframe.grid(row=0, column=0)
        self.rightframe.grid(row=0, column=1)
        self.canvas = FigureCanvasTkAgg(self.fig2, master=self.rightframe)
        self.canvas.get_tk_widget().grid(row=0, column=0)
        
    def fill_leftframe(self):
        pass
    
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
    
    
class Make_Popup_GUI():
    
    def __init__(self, GUI):
        self.GUI = GUI
        self.make_popup()
        self.fill_leftframe()
        self.fill_rightframe()
        
    def make_popup(self):
        self.popup = Toplevel()
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
        
        if not hasattr(self, 'var'):
            # Will already have these attributes if it's being reinitialized
            self.var       = StringVar(value='Input')
            
        Label(frame, text='    Reference Electrode Converter    \n    Coming Soon!    \n').grid(row=0, column=0, columnspan=1)
        Label(frame, text='Variable: ').grid(row=1, column=0, sticky=(E))
        Entry(frame, textvariable=self.var, width=10).grid(row=1, column=1,
                                                                sticky=(W,E))
        return
    
    def fill_rightframe(self):
        # Put relevant buttons in self.leftframe
        frame = self.rightframe
        
        if not hasattr(self, 'res'):
            # Will already have these attributes if it's being reinitialized
            self.res       = StringVar(value='Result')
            
        Label(frame, text='\n\n').grid(row=0, column=0, columnspan=1)
        Label(frame, text='Result: ').grid(row=1, column=0, sticky=(E))
        Entry(frame, textvariable=self.res, width=10).grid(row=1, column=1,
                                                                sticky=(W,E))
        return

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
        
class Popup_Generator():
    '''
    Class to run reference electrode converter
    '''
    
    def __init__(self):
        self.ReferenceElecConvPopup = None
        self.FTdataPopup = None
    
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