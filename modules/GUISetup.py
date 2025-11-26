import tkinter as tk
from tkinter import *
from tkinter.ttk import *
import os
import sys
import numpy as np
import matplotlib.pyplot as plt
from matplotlib.backends.backend_tkagg import FigureCanvasTkAgg

# Get the directory of the current script
script_dir = os.path.dirname(os.path.abspath(__file__))

# Add a directory relative to the script's location
sys.path.append(os.path.join(script_dir, ))
# print(sys.path)

# Now you can import modules
# from modules import Plotter

Overlay_options = ['True', 'False']

bool_map = {"true": True, 'True': True,
            "false": False, 'False': False}

Time_Units = ['hrs','min','s','ms','\u03BCs']

Current_Units = ['A', 'mA', '\u03BCA', 'nA', 'pA', 'fA']

Potential_Units = ['kV', 'V', 'mV', '\u03BCV']

Reference_ = ['SCE', 'Ag/AgCl', 'RHE', 'SHE', 'Ag wire', 'PdH\u2082', 'Fc/Fc\u207A']

Export_format_types = ['png', 'pdf', 'svg', 'jpg']

fig2Options = ['V vs t','I vs t','I vs V',]

EIS_options = ['Nyquist', '|Z| Bode', 'Phase Bode']

Impedance_units = ['G\u03A9', 'M\u03A9', 'k\u03A9', '\u03A9']

Axis_scales = ['linear', 'log', 'logit', 'symlog']

Location_cbar = ['right', 'top']

Legend_loc = ['best', 'upper left', 'upper right',
              'lower left', 'lower right','upper center',
              'lower center', 'center left', 'center right', 'center']

Inset_loc = ['center', 'center left', 'center right',
             'upper center','upper left', 'upper right',
             'lower center','lower left', 'lower right']

mask_options = ['t', 'V', 'I']

# Data for the cascading OptionMenus
cmaps_data = {
    'Sequential 1': ['viridis', 'plasma', 'inferno', 'magma', 'cividis'],
    'Sequential 2': ['Greys', 'Purples', 'Blues', 'Greens', 'Oranges', 'Reds',
                      'YlOrBr', 'YlOrRd', 'OrRd', 'PuRd', 'RdPu', 'BuPu',
                      'GnBu', 'PuBu', 'YlGnBu', 'PuBuGn', 'BuGn', 'YlGn'],
    'Sequential 3': ['binary', 'gist_yarg', 'gist_gray', 'gray', 'bone',
                      'pink', 'spring', 'summer', 'autumn', 'winter', 'cool',
                      'Wistia', 'hot', 'afmhot', 'gist_heat', 'copper'],
    'Diverging': ['PiYG', 'PRGn', 'BrBG', 'PuOr', 'RdGy', 'RdBu', 'RdYlBu',
                      'RdYlGn', 'Spectral', 'coolwarm', 'bwr', 'seismic',
                      'berlin', 'managua', 'vanimo'],
    'Cyclic': ['twilight', 'twilight_shifted', 'hsv'],
    'Qualitative': ['Pastel1', 'Pastel2', 'Paired', 'Accent', 'Dark2',
                      'Set1', 'Set2', 'Set3', 'tab10', 'tab20', 'tab20b',
                      'tab20c'],
    'Miscellaneous': ['flag', 'prism', 'ocean', 'gist_earth', 'terrain',
                      'gist_stern', 'gnuplot', 'gnuplot2', 'CMRmap',
                      'cubehelix', 'brg', 'gist_rainbow', 'rainbow', 'jet',
                      'turbo', 'nipy_spectral', 'gist_ncar',],
    'StyleFile': ['Default']
}
cmaps = ['viridis', 'hot', 'gist_gray', 'afmhot', 'plasma', 'inferno', 
         'magma', 'cividis','Greys', 'Purples', 'Blues', 'Greens', 
         'Oranges', 'Reds', 'YlOrBr', 'YlOrRd', 'OrRd', 'PuRd', 'RdPu', 
         'BuPu','GnBu', 'PuBu', 'YlGnBu', 'PuBuGn', 'BuGn', 'YlGn',
         'PiYG', 'PRGn', 'BrBG', 'PuOr', 'RdGy', 'RdBu', 'RdYlBu',
         'RdYlGn', 'coolwarm', 'bwr']

linestyle_str = ['solid', 'None', 'dotted', 'dashed', 'dashdot']

marker_str = ["None", ".", ",", "o", "v", "^", "<", ">", "1", "2", "3", "4", "8", "s",
              "p", "P", "*", "h", "H", "+", "x", "X", "D", "d", "|", "_", 0, 1, 2, 3]

peak_analysis_options = ['None', 'Reversible', 'Peak Finder']

def where(l, val):
    for i, value in enumerate(l):
        if value == val:
            return i
    raise IndexError


def focus_next_widget(event):
    widget = event.widget.tk_focusNext()
    widget.focus()
    try:
        widget.select_range(0, 'end')
    except:
        pass
    return("break")


def OptionMenuStringVar(frame, options,row,col,sticky,idx=0, return_widget=False):
    var = StringVar()
    menu = OptionMenu(frame, var, options[idx],*options)
    menu.grid(row=row, column=col, sticky=sticky)
    if return_widget:
        return var, menu
    return var


def OptionMenuIntVar(frame, options,row,col,sticky,idx=0, return_widget=False):
    var = IntVar()
    menu = OptionMenu(frame, var, options[idx],*options)
    menu.grid(row=row, column=col, sticky=sticky)
    if return_widget:
        return var, menu
    return var


def EntryStringVar(frame, width, row, col, sticky, default='', bind_key=None,
                   bind_func=None, tab=False, returnTab = False, return_widget=False):
    var = StringVar(value=str(default))
    entry = Entry(frame, width=width, textvariable=var)
    entry.grid(row=row, column=col, sticky=sticky)
    if tab:
        entry.bind('<Tab>', focus_next_widget)
    if returnTab:
        entry.bind('<Return>', focus_next_widget)
    if bind_key:
        entry.bind(bind_key, bind_func)
    if return_widget:
        return var, entry
    return var

def Labels_in_column(frame, labels, column, start_row, sticky):
    widgets = [Label(frame, text=string) for string in labels]
    widgets_in_column(widgets, column, start_row, sticky)
    
        
def widgets_in_column(widgets, column, start_row, sticky):
    row = start_row
    for widget in widgets:
        widget.grid(column=column, row=row, sticky=sticky)
        row += 1


class GUISetupMethods():
    
    def MakeUpdateFrame(self, frame):
        # Reset ADC view button
        Button(frame, text='Update Plot', width=36, 
               command=self.EchemFig.update_plot).grid(
                   row=0, column=0, sticky=(W,E))
        # Copy to Clipboard Button
        Button(frame, text="Copy Figure to Clipboard", width=36,
               command=self.copy_figure_to_clipboard).grid(
                   row=0, column=1, sticky="we", pady=10)
        return
    
    def MakePlotParamsFrame(self, frame):
        Label(frame, text='   Overlay: ').grid(row=0, column=0, sticky=(E))
        self.Overlay_ = OptionMenuStringVar(frame, Overlay_options, 0, 1, (E), idx=1,)
        Label(frame, text='   ').grid(row=0, column=2, sticky=(E))
        Label(frame, text='   File selection: ').grid(row=0, column=3)
        # PointsList selection options
        self.fig2ptselection, self.fig2ptoptmenu = OptionMenuIntVar(frame, [1,], 0, 4, (W),
                                                                    return_widget=True)
        
        Label(frame, text='   Plot type: ').grid(row=1, column=0)
        # Voltammetry view options
        self.fig2selection, self.fig2typeoptmenu = OptionMenuStringVar(
                                                frame, fig2Options,idx=2,
                                                row=1, col=1, sticky=(E),
                                                return_widget=True)
        Label(frame, text='   ').grid(row=1, column=2, sticky=(E))
        Label(frame, text='   EIS type: ').grid(row=1, column=3, sticky=(E), pady=10)
        # EIS view options
        self.EIS_view_selection = OptionMenuStringVar(frame, EIS_options, 1, 4, (E))  
                               
    def MakeEchemFrame(self, frame, ResizeFrame):   
        self.fig = plt.Figure(figsize=(4,4), dpi=150)
        self.fig.add_subplot(111)
        
        self.MakeResizeFrame(ResizeFrame)     
        
        # Canvas setup               
        self.canvas_agg  = FigureCanvasTkAgg(self.fig, master=frame)
        self.canvas_widget = self.canvas_agg.get_tk_widget()
        self.canvas_widget.grid(row=0, column=0, sticky="nsew")
        
        # Allow figure area to expand with window
        frame.rowconfigure(5, weight=1)
        for col in range(4):
            frame.columnconfigure(col, weight=1)
        
        # Bind entry changes
        self.plot_width.trace_add("write", self.on_entry_resize)
        self.plot_height.trace_add("write", self.on_entry_resize)
    
        # Bind window/frame size changes
        self.canvas_widget.bind("<Configure>", self.on_canvas_resize)
        
        # --- FIXED: Allow canvas widget to expand ---
        self.canvas_widget.grid_rowconfigure(0, weight=1)
        self.canvas_widget.grid_columnconfigure(0, weight=1)
    
    def MakeResizeFrame(self, frame):
        # Variables for entry boxes
        self.plot_width = tk.StringVar(value="4")
        self.plot_height = tk.StringVar(value="4")
    
        # Trace callbacks
        self.plot_width.trace_add("write", self.on_entry_resize)
        self.plot_height.trace_add("write", self.on_entry_resize)
    
        # Labels and entries
        Label(frame, text="  Width (in): ").grid(row=3, column=0)
        Entry(frame, textvariable=self.plot_width, width=10).grid(row=3, column=1, sticky=E)
        Label(frame, text='     ').grid(row=3, column=2, sticky=(E))
        Label(frame, text="  Height (in): ").grid(row=3, column=3, pady=10)
        Entry(frame, textvariable=self.plot_height, width=10).grid(row=3, column=4, sticky=E)
    
    def on_canvas_resize(self, event):
        """Resize matplotlib figure to match canvas size exactly."""
        if event.width < 10 or event.height < 10:
            return
        
        dpi = self.fig.get_dpi()
        w_in = event.width / dpi
        h_in = event.height / dpi
    
        self.fig.set_size_inches(w_in, h_in, forward=False)
        self.canvas_agg.draw_idle()
    
    def on_entry_resize(self, *args):
        """User typed a manual width/height → update figure size."""
        try:
            w = float(self.plot_width.get())
            h = float(self.plot_height.get())
            self.fig.set_size_inches(w, h, forward=True)
            self.canvas_agg.draw_idle()
        except ValueError:
            pass
    
    def MakePlotTypeFrame(self, frame):
        tabs = Notebook(frame)
        VoltFrame = Frame(tabs)
        EISFrame = Frame(tabs)
        CustomizeFrame  = Frame(tabs)
        DataControlFrame = Frame(tabs)
        ExportFrame = Frame(tabs)
        
        tabs.add(VoltFrame, text='  Voltamperometric  ')
        tabs.add(EISFrame, text='  EIS  ')
        tabs.add(CustomizeFrame, text='  Customize Plot  ')
        tabs.add(DataControlFrame, text='  Analysis  ')
        tabs.add(ExportFrame, text='  Export  ')
        tabs.pack(expand=1, fill='both')
        
        self.MakeVoltFrame(VoltFrame)
        self.MakeEISFrame(EISFrame)
        self.MakeCustomizeFrame(CustomizeFrame)
        self.MakeDataControlFrame(DataControlFrame)
        self.MakeExportSettingsFrame(ExportFrame)
        
        
    def MakeVoltFrame(self, frame):
        '''Make Frame for Voltamperometric Tab'''
        
        NormFrame = Frame(frame)
        NormFrame.grid(row=0, column=0, sticky=(N,S,E,W), pady=10)
        
        '''Norm'''
        Label(NormFrame, text='Title: ').grid(row=0, column=0, sticky=(E))
        self.title_name = EntryStringVar(NormFrame, 10, 0, 1, (W), tab=True,
                                                  bind_key='<Return>', default='')
        
        Label(NormFrame, text='Time Units: ').grid(row=1, column=0, sticky=(E))
        self.TimeUnits_ = OptionMenuStringVar(NormFrame, Time_Units, 1, 1, (W), idx=2)
        
        Label(NormFrame, text='Reference Type: ').grid(row=2, column=0, sticky=(E))
        self.Reference_update = OptionMenuStringVar(NormFrame, Reference_, 2, 1, (W))
        
        Label(NormFrame, text='Potential Units: ').grid(row=3, column=0, sticky=(E))
        self.PotentialUnits_ = OptionMenuStringVar(NormFrame, Potential_Units, 3, 1, (W), idx=1)
        
        Label(NormFrame, text='Current Units: ').grid(row=4, column=0, sticky=(E))
        self.CurrentUnits_ = OptionMenuStringVar(NormFrame, Current_Units, 4, 1, (W))
        
        Label(NormFrame, text='').grid(row=5, column=0, sticky=(E))
        
        # Make Nested Notebooks
        inner_tabs = Notebook(frame)
        inner_tabs.grid(row=1, column=0, sticky=(N,S,E,W), pady=10)
        frame.grid_rowconfigure(0, weight=1)
        frame.grid_columnconfigure(0, weight=1)
        AxisFrame = Frame(inner_tabs)
        CyclesFrame = Frame(inner_tabs)
        ShiftsFrame = Frame(inner_tabs)
        ColorbarFrame = Frame(inner_tabs)
        
        inner_tabs.add(AxisFrame, text='  Axis  ')
        inner_tabs.add(CyclesFrame, text='  Cycles & Files to Plot  ')
        inner_tabs.add(ShiftsFrame, text='  Shifts  ')
        inner_tabs.add(ColorbarFrame, text='  Colorbar  ')
        
        '''Axis'''
        Label(AxisFrame, text='Box aspect: ').grid(row=0, column=0, sticky=(E))
        self.box_aspect = EntryStringVar(AxisFrame, 10, 0, 1, (W,E), tab=True,
                                                  bind_key='<Return>', default='1')
        Label(AxisFrame, text='').grid(row=1, column=0, sticky=(E))
        Label(AxisFrame, text='X-axis tick multiple: ').grid(row=2, column=0, sticky=(E))
        self.x_axis_tickmultiple = EntryStringVar(AxisFrame, 10, 2, 1, (W,E), tab=True,
                                                  bind_key='<Return>', default='')
        Label(AxisFrame, text='Minor tick multiple: ').grid(row=3, column=0, sticky=(E))
        self.x_axis_minor_tickmultiple = EntryStringVar(AxisFrame, 10, 3, 1, (W,E), tab=True,
                                                  bind_key='<Return>', default='')
        Label(AxisFrame, text='Min: ').grid(row=4, column=0, sticky=(E))
        self.echem_xminval = EntryStringVar(AxisFrame, 10, 4, 1, (W,E), tab=True,
                                                  bind_key='<Return>', default='')
        Label(AxisFrame, text='Max: ').grid(row=4, column=3, sticky=(E))
        self.echem_xmaxval = EntryStringVar(AxisFrame, 10, 4, 4, (W,E), tab=True,
                                                  bind_key='<Return>', default='')
        Label(AxisFrame, text='').grid(row=5, column=0, sticky=(E))
        Label(AxisFrame, text='Y-axis tick multiple: ').grid(row=9, column=0, sticky=(E))
        self.y_axis_tickmultiple = EntryStringVar(AxisFrame, 10, 9, 1, (W,E), tab=True,
                                                  bind_key='<Return>', default='')
        Label(AxisFrame, text='Minor tick multiple: ').grid(row=10, column=0, sticky=(E))
        self.y_axis_minor_tickmultiple = EntryStringVar(AxisFrame, 10, 10, 1, (W,E), tab=True,
                                                  bind_key='<Return>', default='')
        Label(AxisFrame, text='Min: ').grid(row=11, column=0, sticky=(E))
        self.echem_yminval = EntryStringVar(AxisFrame, 10, 11, 1, (W,E), tab=True,
                                                  bind_key='<Return>', default='')
        Label(AxisFrame, text='Max: ').grid(row=11, column=3, sticky=(E))
        self.echem_ymaxval = EntryStringVar(AxisFrame, 10, 11, 4, (W,E), tab=True,
                                                  bind_key='<Return>', default='')
        Label(AxisFrame, text='').grid(row=13, column=0, sticky=(E))
        
        '''Cycles and Files to Plot'''
        Label(CyclesFrame, text='Cycles to Plot: ').grid(row=0, column=0, sticky=(E))
        self.cycles_to_plot = EntryStringVar(CyclesFrame, 10, 0, 1, (W,E), tab=True,
                                                  bind_key='<Return>', default='')
        Label(CyclesFrame, text='Manual Input: ').grid(row=1, column=0, sticky=(E))
        self.cycles_to_plot_manual = OptionMenuStringVar(CyclesFrame, Overlay_options, 1, 1, (W,E), idx=1,)
        Label(CyclesFrame, text='').grid(row=2, column=0, sticky=(E))
        
        Label(CyclesFrame, text='Files to Plot: ').grid(row=3, column=0, sticky=(E))
        self.files_to_plot = EntryStringVar(CyclesFrame, 10, 3, 1, (W,E), tab=True,
                                                  bind_key='<Return>', default='')
        Label(CyclesFrame, text='(Only for Overlay)').grid(row=3, column=2, sticky=(E))
        
        '''Shifts'''
        Label(ShiftsFrame, text='X-axis shifts: ').grid(row=3, column=0, sticky=(E))
        self.x_axis_shifts = EntryStringVar(ShiftsFrame, 10, 3, 1, (W,E), tab=True,
                                                  bind_key='<Return>', default='')
        Label(ShiftsFrame, text='Y-axis shifts: ').grid(row=3, column=2, sticky=(E))
        self.y_axis_shifts = EntryStringVar(ShiftsFrame, 10, 3, 3, (W,E), tab=True,
                                                  bind_key='<Return>', default='')
        Label(ShiftsFrame, text='').grid(row=4, column=0, sticky=(E))
        
        Label(ShiftsFrame, text='Partial Import of data').grid(row=5, column=0, columnspan=2, sticky=(E))
        Label(ShiftsFrame, text='Greater than: ').grid(row=6, column=0, sticky=(E))
        self.start_after_option = OptionMenuStringVar(ShiftsFrame, Overlay_options, 6, 1, (W,E), idx=1,)
        self.start_after_var = OptionMenuStringVar(ShiftsFrame, mask_options, 7, 0, (E), idx=0)
        self.start_after_float = EntryStringVar(ShiftsFrame, 10, 7, 1, (W,E), tab=True,
                                                  bind_key='<Return>', default='')
        
        Label(ShiftsFrame, text='Less than: ').grid(row=8, column=0, sticky=(E))
        self.end_before_option = OptionMenuStringVar(ShiftsFrame, Overlay_options, 8, 1, (W,E), idx=1,)
        self.end_before_var = OptionMenuStringVar(ShiftsFrame, mask_options, 9, 0, (E), idx=0)
        self.end_before_float = EntryStringVar(ShiftsFrame, 10, 9, 1, (W,E), tab=True,
                                                  bind_key='<Return>', default='')
        
        '''Colorbar'''
        Label(ColorbarFrame, text='Colorbar: ').grid(row=0, column=0, sticky=(E))
        self.Colorbar_ = OptionMenuStringVar(ColorbarFrame, Overlay_options, 0, 1, (W,E), idx=1,)
        Label(ColorbarFrame, text='Size: ').grid(row=1, column=0, sticky=(E))
        self.fraction_for_cbar = EntryStringVar(ColorbarFrame, 10, 1, 1, (W,E), tab=True,
                                                  bind_key='<Return>', default='5')
        # Label(ColorbarFrame, text='Shrink: ').grid(row=2, column=0, sticky=(E))
        # self.shrink_for_cbar = EntryStringVar(ColorbarFrame, 10, 2, 1, (W,E), tab=True,
        #                                           bind_key='<Return>', default='1')
        # Label(ColorbarFrame, text='Aspect: ').grid(row=3, column=0, sticky=(E))
        # self.aspect_for_cbar = EntryStringVar(ColorbarFrame, 10, 3, 1, (W,E), tab=True,
        #                                           bind_key='<Return>', default='20')
        Label(ColorbarFrame, text='Pad: ').grid(row=4, column=0, sticky=(E))
        self.pad_for_cbar = EntryStringVar(ColorbarFrame, 10, 4, 1, (W,E), tab=True,
                                                  bind_key='<Return>', default='5')
        
        # Label(ColorbarFrame, text='Label: ').grid(row=15, column=0, sticky=(E))
        # self.cbar_label = EntryStringVar(ColorbarFrame, 10, 15, 1, (W,E), tab=True,
        #                                           bind_key='<Return>', default=' ')
        
        Label(ColorbarFrame, text='Location: ').grid(row=16, column=0, sticky=(E))
        self.Location_for_cbar = OptionMenuStringVar(ColorbarFrame, Location_cbar, 16, 1, (W,E), idx=0,)
        
        Label(ColorbarFrame, text='Tick Labels: ').grid(row=17, column=0, sticky=(E))
        self.labels_for_cbar = EntryStringVar(ColorbarFrame, 10, 17, 1, (W,E), tab=True,
                                                  bind_key='<Return>', default='')
        Label(ColorbarFrame, text='Number of Ticks: ').grid(row=18, column=0, sticky=(E))
        self.ticknum_for_cbar = EntryStringVar(ColorbarFrame, 10, 18, 1, (W,E), tab=True,
                                                  bind_key='<Return>', default='')
        Label(ColorbarFrame, text='').grid(row=19, column=0, sticky=(E))
        
        
        
    def MakeEISFrame(self, frame):
        '''Make Frame for EIS tab'''
        
        NormFrame = Frame(frame)
        NormFrame.grid(row=0, column=0, sticky=(N,S,E,W), pady=10)
        
        '''Norm'''
        Label(NormFrame, text='Title: ').grid(row=0, column=0, sticky=(E))
        self.title_name_EIS = EntryStringVar(NormFrame, 10, 0, 1, (W,E), tab=True,
                                                  bind_key='<Return>', default='')
        Label(NormFrame, text='Impedance Units: ').grid(row=2, column=0, sticky=(E))
        self.ImpedanceUnits_ = OptionMenuStringVar(NormFrame, Impedance_units, 2, 1, (W,E), idx=3)
        Label(NormFrame, text='Scale: ').grid(row=3, column=0, sticky=(E))
        self.Scale_EIS = OptionMenuStringVar(NormFrame, Axis_scales, 3, 1, (W,E), idx=0,)
        
        # Make Nested Notebooks
        inner_tabs = Notebook(frame)
        inner_tabs.grid(row=1, column=0, sticky=(N,S,E,W), pady=10)
        frame.grid_rowconfigure(0, weight=1)
        frame.grid_columnconfigure(0, weight=1)
        AxisFrame = Frame(inner_tabs)
        CyclesFrame = Frame(inner_tabs)
        ShiftsFrame = Frame(inner_tabs)
        
        inner_tabs.add(AxisFrame, text='  Axis  ')
        inner_tabs.add(CyclesFrame, text='  Cycles & Files to Plot  ')
        inner_tabs.add(ShiftsFrame, text='  Shifts  ')
        
        '''Axis'''
        Label(AxisFrame, text='Box aspect: ').grid(row=0, column=0, sticky=(E))
        self.box_aspect_EIS = EntryStringVar(AxisFrame, 10, 0, 1, (W,E), tab=True,
                                                  bind_key='<Return>', default='1')
        Label(AxisFrame, text='').grid(row=1, column=0, sticky=(E))
        Label(AxisFrame, text='X-axis tick multiple: ').grid(row=2, column=0, sticky=(E))
        self.x_axis_tickmultiple_EIS = EntryStringVar(AxisFrame, 10, 2, 1, (W,E), tab=True,
                                                  bind_key='<Return>', default='')
        Label(AxisFrame, text='Minor tick multiple: ').grid(row=3, column=0, sticky=(E))
        self.x_axis_minor_tickmultiple_EIS = EntryStringVar(AxisFrame, 10, 3, 1, (W,E), tab=True,
                                                  bind_key='<Return>', default='')
        Label(AxisFrame, text='Min: ').grid(row=4, column=0, sticky=(E))
        self.echem_xminval_EIS = EntryStringVar(AxisFrame, 10, 4, 1, (W,E), tab=True,
                                                  bind_key='<Return>', default='')
        Label(AxisFrame, text='Max: ').grid(row=4, column=3, sticky=(E))
        self.echem_xmaxval_EIS = EntryStringVar(AxisFrame, 10, 4, 4, (W,E), tab=True,
                                                  bind_key='<Return>', default='')
        Label(AxisFrame, text='').grid(row=5, column=0, sticky=(E))
        Label(AxisFrame, text='Y-axis tick multiple: ').grid(row=10, column=0, sticky=(E))
        self.y_axis_tickmultiple_EIS = EntryStringVar(AxisFrame, 10, 10, 1, (W,E), tab=True,
                                                  bind_key='<Return>', default='')
        Label(AxisFrame, text='Minor tick multiple: ').grid(row=11, column=0, sticky=(E))
        self.y_axis_minor_tickmultiple_EIS = EntryStringVar(AxisFrame, 10, 11, 1, (W,E), tab=True,
                                                  bind_key='<Return>', default='')
        Label(AxisFrame, text='Min: ').grid(row=12, column=0, sticky=(E))
        self.echem_yminval_EIS = EntryStringVar(AxisFrame, 10, 12, 1, (W,E), tab=True,
                                                  bind_key='<Return>', default='')
        Label(AxisFrame, text='Max: ').grid(row=12, column=3, sticky=(E))
        self.echem_ymaxval_EIS = EntryStringVar(AxisFrame, 10, 12, 4, (W,E), tab=True,
                                                  bind_key='<Return>', default='')
        Label(frame, text='').grid(row=13, column=0, sticky=(E))
        
        '''Cycles and Files to Plot'''
        Label(CyclesFrame, text='Cycles to Plot: ').grid(row=0, column=0, sticky=(E))
        self.cycles_to_plot_EIS = EntryStringVar(CyclesFrame, 10, 0, 1, (W,E), tab=True,
                                                  bind_key='<Return>', default='')
        Label(CyclesFrame, text='Manual Input: ').grid(row=1, column=0, sticky=(E))
        self.cycles_to_plot_manual_EIS = OptionMenuStringVar(CyclesFrame, Overlay_options, 1, 1, (W,E), idx=1,)
        Label(CyclesFrame, text='').grid(row=2, column=0, sticky=(E))
        
        Label(CyclesFrame, text='Files to Plot: ').grid(row=3, column=0, sticky=(E))
        self.files_to_plot_EIS = EntryStringVar(CyclesFrame, 10, 3, 1, (W,E), tab=True,
                                                  bind_key='<Return>', default='')
        Label(CyclesFrame, text='(Only for Overlay)').grid(row=3, column=2, sticky=(E))
        
        '''Shifts'''
        Label(ShiftsFrame, text='X-axis shifts: ').grid(row=0, column=0, sticky=(E))
        self.x_axis_shifts_EIS = EntryStringVar(ShiftsFrame, 10, 0, 1, (W,E), tab=True,
                                                  bind_key='<Return>', default='')
        Label(ShiftsFrame, text='Y-axis shifts: ').grid(row=0, column=2, sticky=(E))
        self.y_axis_shifts_EIS = EntryStringVar(ShiftsFrame, 10, 0, 3, (W,E), tab=True,
                                                  bind_key='<Return>', default='')
        
    def MakeCustomizeFrame(self, frame):
        '''Make Frame for Color, Legend, Lines, and Markers'''
        
        ColorFrame = Frame(frame)
        ColorFrame.grid(row=0, column=0, sticky=(N,S,E,W), pady=10)
        
        '''Color'''
        Label(ColorFrame, text='Color: ').grid(row=1, column=0, sticky=(E))
        
        # First OptionMenu (Category)
        self.category_var = StringVar(ColorFrame)
        # self.category_var.set("Select Category")  # Default category
        
        self.category_menu = OptionMenu(ColorFrame, self.category_var,
                                        list(cmaps_data.keys())[0],
                                        *list(cmaps_data.keys()),
                                        command=self.update_items)
        self.category_menu.grid(row=1, column=1, sticky=(E))

        # Second OptionMenu (Items within Category)
        self.plot_cmap = StringVar(ColorFrame)
        self.cmap_options = OptionMenu(ColorFrame, self.plot_cmap,'',
                                       command=self.EchemFig.update_plot)
        self.cmap_options.grid(row=1, column=3, sticky=(E))  # Initial empty OptionMenu

        # Initially populate the items menu based on the default category
        self.update_items(self.category_var.get())
        
        self.echem_cmap_minval = StringVar(value='0')
        self.echem_cmap_maxval = StringVar(value='1')
        
        Label(ColorFrame, text='Min: ').grid(row=2, column=0, sticky=(E))
        self.echem_cmap_minval = EntryStringVar(ColorFrame, 15, 2, 1, (W,E), tab=True,
                                                  bind_key='<Return>', bind_func=
                                                  self.EchemFig.update_colormap(),
                                                  default='0')
        Label(ColorFrame, text='Max: ').grid(row=2, column=3, sticky=(E))
        self.echem_cmap_maxval = EntryStringVar(ColorFrame, 15, 2, 4, (W,E), tab=True,
                                                  bind_key='<Return>', bind_func=
                                                  self.EchemFig.update_colormap(),
                                                  default='1')
        
        # Make Nested Notebooks
        inner_tabs = Notebook(frame)
        inner_tabs.grid(row=1, column=0, sticky=(N,S,E,W), pady=10)
        frame.grid_rowconfigure(0, weight=1)
        frame.grid_columnconfigure(0, weight=1)
        LegendLSMSFrame = Frame(inner_tabs)
        InsetFrame = Frame(inner_tabs)
        
        inner_tabs.add(LegendLSMSFrame, text='  Legend, Linestyle, and Makers  ')
        inner_tabs.add(InsetFrame, text='  Inset  ')
        
        '''Legend'''
        Label(LegendLSMSFrame, text='Legend: ').grid(row=0, column=0, sticky=(E))
        self.Legend_ = OptionMenuStringVar(LegendLSMSFrame, Overlay_options, 0, 1, (W,E), idx=1,)
        
        Label(LegendLSMSFrame, text='Labels: ').grid(row=1, column=0, sticky=(E))
        self.labels_for_legend = EntryStringVar(LegendLSMSFrame, 10, 1, 1, (W,E), tab=True,
                                                  bind_key='<Return>', default='')
        
        Label(LegendLSMSFrame, text='Location: ').grid(row=2, column=0, sticky=(E))
        self.location_for_legend = OptionMenuStringVar(LegendLSMSFrame, Legend_loc, 2, 1, (W,E), idx=0,)
        
        Label(LegendLSMSFrame, text='Size: ').grid(row=3, column=0, sticky=(E))
        self.fontsize_for_legend = EntryStringVar(LegendLSMSFrame, 10, 3, 1, (W,E), tab=True,
                                                  bind_key='<Return>', default='')
        
        Label(LegendLSMSFrame, text='Handle Length: ').grid(row=7, column=0, sticky=(E))
        self.handle_length_legend = EntryStringVar(LegendLSMSFrame, 10, 7, 1, (W,E), tab=True,
                                                  bind_key='<Return>', default='1.5')
        Label(LegendLSMSFrame, text='').grid(row=8, column=0, sticky=(E))
        
        '''Linestyle and Makers'''
        Label(LegendLSMSFrame, text='Line width: ').grid(row=9, column=0, sticky=(E))
        self.line_width = EntryStringVar(LegendLSMSFrame, 10, 9, 1, (W,E), tab=True,
                                                  bind_key='<Return>', default='2')
        
        Label(LegendLSMSFrame, text='Line Style: ').grid(row=10, column=0, sticky=(E))
        self.line_style = OptionMenuStringVar(LegendLSMSFrame, linestyle_str, 10, 1, (W,E), idx=0,)
        
        Label(LegendLSMSFrame, text='Manual LS: ').grid(row=11, column=0, sticky=(E))
        self.manual_line_styles = EntryStringVar(LegendLSMSFrame, 10, 11, 1, (W,E), tab=True,
                                                  bind_key='<Return>', default='')
        
        Label(LegendLSMSFrame, text='').grid(row=12, column=0, sticky=(E))
        
        Label(LegendLSMSFrame, text='Marker size: ').grid(row=13, column=0, sticky=(E))
        self.marker_width = EntryStringVar(LegendLSMSFrame, 10, 13, 1, (W,E), tab=True,
                                                  bind_key='<Return>', default='5')
        
        Label(LegendLSMSFrame, text='Marker Style: ').grid(row=14, column=0, sticky=(E))
        self.marker_style = OptionMenuStringVar(LegendLSMSFrame, marker_str, 14, 1, (W,E), idx=0,)
        
        Label(LegendLSMSFrame, text='Manual MS: ').grid(row=15, column=0, sticky=(E))
        self.manual_marker_styles = EntryStringVar(LegendLSMSFrame, 10, 15, 1, (W,E), tab=True,
                                                  bind_key='<Return>', default='')
        
        '''Inset'''
        Label(InsetFrame, text='Inset: ').grid(row=0, column=0, sticky=(E))
        self.Inset_ = OptionMenuStringVar(InsetFrame, Overlay_options, 0, 1, (W,E), idx=1,)
        
        Label(InsetFrame, text='% Width: ').grid(row=1, column=0, sticky=(E))
        self.percent_width_Inset = EntryStringVar(InsetFrame, 10, 1, 1, (W,E), tab=True,
                                                  bind_key='<Return>', default='25')
        Label(InsetFrame, text='% Height: ').grid(row=1, column=2, sticky=(E))
        self.percent_height_Inset = EntryStringVar(InsetFrame, 10, 1, 3, (W,E), tab=True,
                                                  bind_key='<Return>', default='25')
        
        Label(InsetFrame, text='Location: ').grid(row=2, column=0, sticky=(E))
        self.location_for_Inset = OptionMenuStringVar(InsetFrame, Inset_loc, 2, 1, (W,E), idx=0,)
        
        Label(InsetFrame, text='Color: ').grid(row=3, column=0, sticky=(E))
        self.Color_Inset = EntryStringVar(InsetFrame, 10, 3, 1, (W,E), tab=True,
                                                  bind_key='<Return>', default='red')
        
        Label(InsetFrame, text='').grid(row=4, column=0, sticky=(E))
        Label(InsetFrame, text='Ticks: ').grid(row=5, column=0, sticky=(E))
        self.Ticks_Inset = OptionMenuStringVar(InsetFrame, Overlay_options, 5, 1, (W,E), idx=1,)
        
        Label(InsetFrame, text='X-axis tick multiple: ').grid(row=6, column=0, sticky=(E))
        self.x_axis_tickmultiple_Inset = EntryStringVar(InsetFrame, 10, 6, 1, (W,E), tab=True,
                                                  bind_key='<Return>', default='')
        Label(InsetFrame, text='Min: ').grid(row=7, column=0, sticky=(E))
        self.xminval_Inset = EntryStringVar(InsetFrame, 10, 7, 1, (W,E), tab=True,
                                                  bind_key='<Return>', default='0')
        Label(InsetFrame, text='Max: ').grid(row=7, column=2, sticky=(E))
        self.xmaxval_Inset = EntryStringVar(InsetFrame, 10, 7, 3, (W,E), tab=True,
                                                  bind_key='<Return>', default='1')     
        
        Label(InsetFrame, text='Y-axis tick multiple: ').grid(row=10, column=0, sticky=(E))
        self.y_axis_tickmultiple_Inset = EntryStringVar(InsetFrame, 10, 10, 1, (W,E), tab=True,
                                                  bind_key='<Return>', default='')
        Label(InsetFrame, text='Min: ').grid(row=12, column=0, sticky=(E))
        self.yminval_Inset = EntryStringVar(InsetFrame, 10, 12, 1, (W,E), tab=True,
                                                  bind_key='<Return>', default='0')
        Label(InsetFrame, text='Max: ').grid(row=12, column=2, sticky=(E))
        self.ymaxval_Inset = EntryStringVar(InsetFrame, 10, 12, 3, (W,E), tab=True,
                                                  bind_key='<Return>', default='1')
        Label(InsetFrame, text='').grid(row=13, column=0, sticky=(E))
        
        Label(InsetFrame, text='Spine Line Style: ').grid(row=14, column=0, sticky=(E))
        self.spine_line_style_Inset = OptionMenuStringVar(InsetFrame, linestyle_str, 14, 1, (W,E), idx=0,)
        Label(InsetFrame, text='Zoom Line Style: ').grid(row=15, column=0, sticky=(E))
        self.zoom_line_style_Inset = OptionMenuStringVar(InsetFrame, linestyle_str, 15, 1, (W,E), idx=0,)
        Label(InsetFrame, text='Connecting Line Styles').grid(row=16, column=0, sticky=(E))
        
        Label(InsetFrame, text='Upper Left: ').grid(row=17, column=0, sticky=(E))
        self.upperleft_line_style_Inset = OptionMenuStringVar(InsetFrame, linestyle_str, 17, 1, (W,E), idx=1,)
        Label(InsetFrame, text='Upper Right: ').grid(row=17, column=2, sticky=(E))
        self.upperright_line_style_Inset = OptionMenuStringVar(InsetFrame, linestyle_str, 17, 3, (W,E), idx=1,)
        Label(InsetFrame, text='Lower Left: ').grid(row=18, column=0, sticky=(E))
        self.lowerleft_line_style_Inset = OptionMenuStringVar(InsetFrame, linestyle_str, 18, 1, (W,E), idx=1,)
        Label(InsetFrame, text='Lower Right: ').grid(row=18, column=2, sticky=(E))
        self.lowerright_line_style_Inset = OptionMenuStringVar(InsetFrame, linestyle_str, 18, 3, (W,E), idx=1,)
        
    def MakeDataControlFrame(self, frame):
        Label(frame, text='Apply Notch Filter: ').grid(row=0, column=0, sticky=(E))
        self.apply_notch_filter = OptionMenuStringVar(frame, Overlay_options, 0, 1, (W,E), idx=1,)
        
        Label(frame, text='Frequencies: ').grid(row=1, column=0, sticky=(E))
        self.freqs_for_notch_filter = EntryStringVar(frame, 10, 1, 1, (W,E), tab=True,
                                                  bind_key='<Return>', default='')
        Label(frame, text='').grid(row=2, column=0, sticky=(E))
        Label(frame, text='Peak Analysis: ').grid(row=6, column=0, sticky=(E))
        self.AnalyzePeak_ = OptionMenuStringVar(frame, peak_analysis_options, 6, 1, (W,E), idx=0,)
        Label(frame, text='Prominence: ').grid(row=7, column=0, sticky=(E))
        self.prominence_peak = EntryStringVar(frame, 10, 7, 1, (W), tab=True,
                                                  bind_key='<Return>', default='1e-2')
        Label(frame, text='Height: ').grid(row=7, column=2, sticky=(E))
        self.height_peak = EntryStringVar(frame, 10, 7, 3, (W), tab=True,
                                                  bind_key='<Return>', default='1e-2')
        Label(frame, text='Oxidative Volt: ').grid(row=8, column=0, sticky=(E))
        self.oxidative_baseline_voltage = EntryStringVar(frame, 10, 8, 1, (W), tab=True,
                                                  bind_key='<Return>', default='')
        Label(frame, text='Reductive Volt: ').grid(row=8, column=2, sticky=(E))
        self.reductive_baseline_voltage = EntryStringVar(frame, 10, 8, 3, (W), tab=True,
                                                  bind_key='<Return>', default='')
        Label(frame, text='Peak Color: ').grid(row=9, column=0, sticky=(E))
        self.peak_analysis_color = EntryStringVar(frame, 10, 9, 1, (W), tab=True,
                                                  bind_key='<Return>', default='red')
        Label(frame, text='Baseline Color: ').grid(row=9, column=2, sticky=(E))
        self.baseline_analysis_color = EntryStringVar(frame, 10, 9, 3, (W), tab=True,
                                                  bind_key='<Return>', default='green')
        Label(frame, text='Line width: ').grid(row=10, column=0, sticky=(E))
        self.baseline_analysis_line_width = EntryStringVar(frame, 10, 10, 1, (W), tab=True,
                                                  bind_key='<Return>', default='1')
        Label(frame, text='Line Style: ').grid(row=10, column=2, sticky=(E))
        self.baseline_analysis_line_style = OptionMenuStringVar(frame, linestyle_str, 10, 3, (W,E), idx=3,)
    
    def update_items(self, selected_category):
        # Clear existing options in the item menu
        self.cmap_options['menu'].delete(0, 'end')

        # Get the new items based on the selected category
        new_items = cmaps_data.get(selected_category, [])

        # Add new items to the item menu
        for item in new_items:
            self.cmap_options['menu'].add_command(label=item, command=tk._setit(self.plot_cmap, item))

        # Set the default for the plot_cmap to the first item in the new list, or an empty string if no items
        self.plot_cmap.set(new_items[0] if new_items else "")
        
    def MakeExportSettingsFrame(self, frame):
        #Trans - True/False
        #dpi #
        #format - 'png', 'pdf', 'svg', 'jpg'
        Label(frame, text='Transparent: ').grid(row=1, column=0, sticky=(E))
        self.Transparency = OptionMenuStringVar(frame, Overlay_options, 1, 1, (W,E), idx=1,)
        
        Label(frame, text='dpi: ').grid(row=2, column=0, sticky=(E))
        self.dpi = EntryStringVar(frame, 10, 2, 1, (W,E), tab=True,
                                                  bind_key='<Return>', default='300')
        
        Label(frame, text='Format: ').grid(row=3, column=0, sticky=(E))
        self.Export_format = OptionMenuStringVar(frame, Export_format_types, 3, 1, (W,E), idx=0)
        
        Label(frame, text='').grid(row=4, column=0, sticky=(E))