import sys
import os
from matplotlib import pyplot as plt, cm, ticker as tk
import numpy as np
from tkinter import *
from tkinter.ttk import *
from tkinter import filedialog, messagebox
import time


# Get the directory of the current script
script_dir = os.path.dirname(os.path.abspath(__file__))

# Add a directory relative to the script's location
sys.path.append(os.path.join(script_dir, ))
# print(sys.path)

# Now you can import modules
from modules.Extract_data import extract_data
from modules.GUISetup import GUISetupMethods
from modules.Plotter import EchemFig
from modules.Extras import Popup_Generator

default_stdout = sys.stdout
default_stdin  = sys.stdin
default_stderr = sys.stderr

bool_map = {"true": True, 'True': True,
            "false": False, 'False': False}

class PrintLogger(): 
    '''
    File like object to print console output into Tkinter window
    set sys.stdout = PrintLogger, then print() will print to
    PrintLogger.textbox
    '''
    def __init__(self, textbox): 
        self.textbox = textbox # tk.Text object
        self.textbox.tag_config("red", foreground="red")
        self.textbox.tag_config('black', foreground='black')
        self.textbox.tag_config('blue', foreground='blue')
        self.textbox.tag_config('green', foreground='green')

    def write(self, text):
        try:
            self.textbox.configure(state='normal')
            color='black'
            for msg in ('Warning', 'warning', 'Error', 'error','Input', 'Invalid'):
                if msg in text:
                    color='red'
            if 'Ploting' in text:
                color='blue'
            if 'exported' in text:
                color='green'
            self.textbox.insert('end', text, color) # write text to textbox
            self.textbox.see('end') # scroll to end
            self.textbox.configure(state='disabled')
        except (TclError, RuntimeError):
            # Widget has likely been destroyed — ignore the print
            pass

    def flush(self): # needed for file like object
        pass
    
class GUI(GUISetupMethods, EchemFig, extract_data):
    '''
    Graphical user interface
    '''
    def __init__(self):
        self.root = root
        self.params = {} # master dict to store all parameters
                
        root.title("Plot echem!")
        root.attributes('-topmost', 1)
        root.attributes('-topmost', 0) 
        root.option_add('*tearOff', False)
        
        
        # Menu bar
        menubar         = Menu(root)
        root['menu']    = menubar
        
        menu_file       = Menu(menubar)
        menubar.add_cascade(menu=menu_file, label='File')
        menu_file.add_command(label='Open folder', command=self.openFolder)
        menu_file.add_command(label='Open one file...', command=self.openFile)
        menu_file.add_command(label='Open MPL Style', command=self.openStyleFile)
        menu_file.add_command(label='Export Plot', command=self.save)
        
        menu_fun = Menu(menubar)
        menubar.add_cascade(menu=menu_fun, label='Extra Tools')
        menu_fun.add_command(label='Open Reference Electrode Converter', command=self.ref_converter)
        menu_fun.add_command(label='Fourier Transform Data', command=self.FT_data)
        
        ### SET UP FRAMES ###
        
        # Root has 2 rows: top content + console
        self.root.rowconfigure(0, weight=1)
        self.root.rowconfigure(1, weight=0)
        self.root.columnconfigure(0, weight=0)
        self.root.columnconfigure(1, weight=1)
        
        leftpanel  = Frame(self.root)
        rightpanel = Frame(self.root)
        ConsoleFrame = Frame(self.root)
        
        leftpanel.grid(row=0, column=0, sticky=(N,S,E,W))
        rightpanel.grid(row=0, column=1, sticky=(N,S,E,W))
        ConsoleFrame.grid(row=1, column=0, columnspan=2, sticky=(N,S,E,W))
        
        # Inside console frame
        console = Text(ConsoleFrame, width=97, height=10)
        console.grid(row=0, column=0, sticky=(N,S,E,W))
        
        # --- Input line ---
        self.entry = Entry(ConsoleFrame, width=97)
        self.entry.grid(row=1, column=0, sticky=(N,S,E,W), padx=5, pady=5)
        self.entry.bind("<Return>", self._on_enter)

        # Variable to capture "input()" style responses
        self._input_var = None
        self._input_value = None
        
        # Inside left panel
        PlotTypeFrame = Frame(leftpanel)
        SettingsFrame   = Frame(leftpanel)
        
        PlotTypeFrame.grid(row=2, column=0, sticky=(N,S,E,W))
        SettingsFrame.grid(row=3, column=0, sticky=(N,S,E,W))
        
        # Inside right panel — use grid consistently here
        rightpanel.rowconfigure(0, weight=1)
        rightpanel.columnconfigure(1, weight=1)
        
        UpdateButtonFrame = Frame(rightpanel)
        PlotParamsFrame = Frame(rightpanel)
        EchemFrame   = Frame(rightpanel)
        ResizeFrame = Frame(rightpanel)
        ZoomFrame = Frame(rightpanel)
        
        UpdateButtonFrame.grid(row=0, column=1, sticky=(N,S,E,W))
        PlotParamsFrame.grid(row=1, column=1, sticky=(N,S,E,W))
        ResizeFrame.grid(row=2, column=1, sticky=(N,S,E,W))
        EchemFrame.grid(row=3, column=1, sticky=(N,S,E,W))
        ZoomFrame.grid(row=4, column=1, sticky=(N,S,E,W))
        
        EchemFrame.rowconfigure(0, weight=1)
        EchemFrame.columnconfigure(0, weight=1)
        
        # All inherited from GUISetupMethods
        self.MakeEchemFrame(EchemFrame, ResizeFrame)
            
        # Initialize plotter
        self.EchemFig = EchemFig(self.fig, self)
        self.MakePlotParamsFrame(PlotParamsFrame)
        self.MakeUpdateFrame(UpdateButtonFrame)
        
        self.MakePlotTypeFrame(PlotTypeFrame)
        self.MakeSettingsFrame(SettingsFrame)
        
        self.fig2selection.trace('w', self.fig_opt_changed)
        self.fig2ptselection.trace('w', self.fig_opt_changed)
        
        # Norm traces
        self.CurrentUnits_.trace('w', self.fig_opt_changed)
        self.Reference_update.trace('w', self.fig_opt_changed)
        self.PotentialUnits_.trace('w', self.fig_opt_changed)
        self.TimeUnits_.trace('w', self.fig_opt_changed)
        
        # Both
        self.Overlay_.trace('w', self.fig_opt_changed)
        self.plot_cmap.trace('w', self.fig_opt_changed)
        self.Colorbar_.trace('w', self.fig_opt_changed)
        
        # EIS traces
        self.EIS_view_selection.trace('w', self.fig_opt_changed)
        self.ImpedanceUnits_.trace('w', self.fig_opt_changed)
        self.Scale_EIS.trace('w', self.fig_opt_changed)
        
        # Legend traces
        self.Legend_.trace('w', self.fig_opt_changed)
        self.location_for_legend.trace('w', self.fig_opt_changed)
        
        # Lines and Marker traces
        self.line_style.trace('w', self.fig_opt_changed)
        self.marker_style.trace('w', self.fig_opt_changed)
        
        # Data control traces
        self.apply_notch_filter.trace('w', self.setting_notch_filter)
        
        # Partical import data traces
        self.start_after_option.trace('w', self.fig_opt_changed)
        self.end_before_option.trace('w', self.fig_opt_changed)
        
        # Send print messages to the console
        sys.stdout = PrintLogger(console)
        
    #################### END __init__ ##############################
    
    ########## GUI CALLBACKS ###########    
    
    def openFolder(self):
        folder = filedialog.askdirectory()
        if folder == '':
            print('Error: Folder not chosen')
            return
        # folder = filedialog.askdirectory(initialdir='Z:\Projects\Miguel\Raw data')
        # print(f'Folder path: {folder}')
        try:
            extracted_data, file_max = self.read_file(folder, True, self)
        except TypeError:
            return
        # print(file_max, extracted_data[0], extracted_data[1])
        if len(extracted_data) == 0:
            print('Error: No file .txt, .csv, .mat, .asc')
            return
        self.EchemFig.initialize()
        self.update_fig_to_plot(file_max)
        # if bool_map.get(self.Overlay_.get().strip().lower(), False) == False:
        #     print('Plotting first file...')
        self.EchemFig.set_data_from_file(extracted_data, file_max)
        self.EchemFig.update_plot()
        return
        
    def openFile(self):
        f = filedialog.askopenfilename()
        if f == '':
            print('Error: File not chosen')
            return

        Multi_files = False
        Overlay = bool_map.get(self.Overlay_.get().strip().lower(), False)
        if Overlay == True:
            print('Cannot Overlay one plot: Set Overlay to "False"')
            return
        
        try:
            extracted_data, file_max = self.read_file(f, Multi_files, self)
        except TypeError:
            return
        # self.EchemFig.initialize()

        self.EchemFig.set_data_from_file([extracted_data], 1)
        self.EchemFig.update_plot()
    
    def openStyleFile(self):
        Style_File_dir = filedialog.askopenfilename()
        if Style_File_dir.endswith('.mplstyle'):
            self.EchemFig.load_style_file(Style_File_dir)
        else:
            print('Error: Style File (.mplstyle) not loaded.')
    
    def _on_enter(self, event=None):
        if self._input_var is not None:
            self._input_value = self.entry.get().strip()
            self.entry.delete(0, 'end')
            self._input_var.set(True)   # releases wait_variable()
            
    def console_input(self, prompt=""):
        """Replacement for input()."""
        print(prompt, end="")  # show prompt in console
        self._input_var = BooleanVar()
        self._input_value = None
        self.entry.focus_set()
        self.root.wait_variable(self._input_var)  # wait until user presses Enter
        return self._input_value
            
    def exit_(self):
        # Shutdown procedure
        sys.stdout = default_stdout
        self.root.destroy()
        # sys.exit(0)
        print("Closing GUI...")
        return
        
    ########## DISPLAY FIGURE CALLBACKS ###########
    
    # Selected new view for fig
    def fig_opt_changed(self, *args):
        self.EchemFig.update_plot()
        return
    
    def update_fig_to_plot(self, file_max):
        '''
        Update Echem figure selection dropdowns based on how many files 
        are in the folder and what type is currently selected
        '''
        # Set dropdown to select DataPoint from PointsList
        menu_length = self.fig2ptoptmenu['menu'].index("end") + 1
        if (file_max > 1) and (file_max != menu_length):
            # print(f'Detected {file_max} files, currently {menu_length} in menu')
            menu = self.fig2ptoptmenu['menu']
            menu.delete(0, 'end')
            opts = [i+1 for i in range(file_max)]
            self.fig2ptoptmenu.set_menu(opts[0], *opts)
    
    ########## Extra Functionality CALLBACKS ###########
    
    def ref_converter(self):
        if not hasattr(self, 'Popup_Generator'):
            self.Popup_Generator = Popup_Generator()
        self.Popup_Generator.get(self, 'Ref_converter', None)
        
    def FT_data(self):
        # Get data from the lines
        x_data, y_data = [], []

        for line in self.EchemFig.ax.get_lines():
            x = line.get_xdata()
            y = line.get_ydata()
            if len(x) > 0 and len(y) > 0:
                x_data.append(x)
                y_data.append(y)
        
        # if len(data['x']) == 0:
        #     print('Error: Load file first')
        # else:
        #     if not hasattr(self, 'Popup_Generator'):
        #         self.Popup_Generator = Popup_Generator()
        #     self.Popup_Generator.get(self, 'FT_data', data)
    
    def setting_notch_filter(self, *args):
        print('Warning: reload data file to apply or remove filter')
        return
                
if __name__ == '__main__':
    root = Tk()
    gui = GUI()
    root.protocol("WM_DELETE_WINDOW", gui.exit_)
    root.mainloop()
