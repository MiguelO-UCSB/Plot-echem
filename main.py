import sys
import os
from matplotlib import pyplot as plt, cm, ticker as tk
import numpy as np
from tkinter import *
from tkinter.ttk import *
from tkinter import filedialog, messagebox, font
import time
import io
import subprocess
from PIL import Image

# Windows needs pywin32
if sys.platform.startswith("win"):
    import win32clipboard
    
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
            if msg in ('Ploting', 'plotting', 'Calculated', 'calculated'):
                if msg in text:
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

        # Creating a Font object of "TkDefaultFont"
        self.defaultFont = font.nametofont("TkDefaultFont")

        # Overriding default-font with custom settings
        # i.e changing font-family, size and weight
        root.option_add("*font", self.defaultFont)
        self.defaultFont.configure(family=self.defaultFont['family'],
                                   size=self.defaultFont['size'],
                                   weight=self.defaultFont['weight'])
        
        self.boldFont = font.Font(family=self.defaultFont['family'],
                                  size=self.defaultFont['size'],
                                  weight='bold')
        self.boldStyle = Style().configure('Bold.TButton', font=(self.defaultFont['family'],
                                                               self.defaultFont['size'],
                                                               'bold'))
        self.notebookStyle = Style().configure('TNotebook.Tab', font=(self.defaultFont['family'],
                                                                self.defaultFont['size'],
                                                                'bold'))
        # ---- ICON LOADING
        if sys.platform.startswith("win"):
            # Windows requires .ico
            logo_path = os.path.join(script_dir, "assets/GUI logo.ico")
            root.iconbitmap(logo_path)
        
        elif sys.platform == "darwin" or sys.platform.startswith("linux"):
            # PNG works for macOS + Linux
            try:
                logo_path = os.path.join(script_dir, "assets/GUI logo.png")
                photo = PhotoImage(file=logo_path)
                root.iconphoto(False, photo)
            except TclError:
                print("icon file not loaded")
        
        else:
            print("⚠ Icon not supported on this platform.")
        
        # ---- MENU BAR
        menubar         = Menu(root)
        root['menu']    = menubar
        
        menu_file       = Menu(menubar)
        menubar.add_cascade(menu=menu_file, label='File')
        menu_file.add_command(label='Open folder', command=self.openFolder)
        menu_file.add_command(label='Open one file...', command=self.openFile)
        menu_file.add_command(label='Open MPL Style', command=self.openStyleFile)
        menu_file.add_command(label='Save Plot', command=self.save)
        
        menu_fun = Menu(menubar)
        menubar.add_cascade(menu=menu_fun, label='Extra Tools')
        menu_fun.add_command(label='Open Reference Electrode and Unit Converter', command=self.ref_converter)
        menu_fun.add_command(label='Open Macrodisk RS Calculator', command=self.macro_RS)
        menu_fun.add_command(label='Open Levich Calculator', command=self.macro_RS_RDE)
        menu_fun.add_command(label='Open Microdisk RS Calculator', command=self.micro_RS)
        menu_fun.add_command(label='Fourier Transform Data', command=self.FT_data)

        # ---- MAIN LAYOUT      
        # Root window stretches
        root.rowconfigure(0, weight=1)
        root.columnconfigure(0, weight=1)
        
        # ======= PANED WINDOW (Horizontal split) =======
        main_pane = PanedWindow(root, orient="horizontal")
        main_pane.grid(row=0, column=0, sticky="nsew")
        
        # ---- LEFT PANEL (Now Scrollable)
        
        # 1. Create a container frame to hold both the Canvas and the Scrollbar
        left_container = Frame(main_pane)
        left_container.rowconfigure(0, weight=1)
        left_container.columnconfigure(0, weight=1)
        
        # 2. Create the Canvas and the Scrollbar
        left_canvas = Canvas(left_container)
        left_canvas.grid(row=0, column=0, sticky="nsew")
        
        left_scrollbar = Scrollbar(left_container, orient="vertical", command=left_canvas.yview)
        left_scrollbar.grid(row=0, column=1, sticky="ns")
        
        left_canvas.configure(yscrollcommand=left_scrollbar.set)
        
        # 3. Create your original leftpanel INSIDE the canvas
        leftpanel = Frame(left_canvas)
        
        # 4. Add the leftpanel to a "window" inside the canvas
        # We save the frame_id so we can force it to match the canvas width later
        frame_id = left_canvas.create_window((0, 0), window=leftpanel, anchor="nw")
        
        # 5. Bind events to update scrolling and width dynamically
        # Update the scroll region whenever the leftpanel changes size
        leftpanel.bind(
            "<Configure>",
            lambda e: left_canvas.configure(scrollregion=left_canvas.bbox("all"))
        )
        
        # Force the leftpanel to stretch to the width of the canvas so it doesn't collapse
        left_canvas.bind(
            "<Configure>",
            lambda e: left_canvas.itemconfig(frame_id, width=max(e.width, 1000))
        )
        
        # --- Your original sub-frames (unchanged, just nested in the new leftpanel) ---
        UpdateButtonFrame = Frame(leftpanel)
        PlotParamsFrame = Frame(leftpanel)
        ResizeFrame = Frame(leftpanel)
        PlotTypeFrame = Frame(leftpanel)
        
        UpdateButtonFrame.grid(row=0, column=0, sticky=(N,S,E,W))
        PlotParamsFrame.grid(row=1, column=0, sticky=(N,S,E,W))
        ResizeFrame.grid(row=2, column=0, sticky=(N,S,E,W))
        PlotTypeFrame.grid(row=3, column=0, sticky=(N,S,E,W))
        
        # Add the CONTAINER to the paned window, not the inner panel
        main_pane.add(left_container)   
        
        # ---- RIGHT PANEL
        rightpanel = Frame(main_pane)
        rightpanel.rowconfigure(0, weight=1)  # EchemFrame grows
        rightpanel.columnconfigure(0, weight=1)
        
        EchemFrame   = Frame(rightpanel)
        EchemFrame.grid(row=0, column=0, sticky=(N,S,E,W))
        
        main_pane.add(rightpanel)
        
        # ---- CONSOLE FRAME
        ConsoleFrame = Frame(root)
        ConsoleFrame.grid(row=1, column=0, sticky="nsew")
        
        ConsoleFrame.columnconfigure(0, weight=1)
        ConsoleFrame.rowconfigure(0, weight=1)
        
        console = Text(ConsoleFrame, width=100, height=10)
        console.grid(row=0, column=0, sticky=(N,S,E,W))
        
        self.entry = Entry(ConsoleFrame)
        self.entry.grid(row=1, column=0, sticky=(N,S,E,W), padx=5, pady=5)
        self.entry.bind("<Return>", self._on_enter)

        # Variable to capture "input()" style responses
        self._input_var = None
        self._input_value = None
        
        # ---- GUI SETUP & TRACES
        # All inherited from GUISetupMethods
        self.MakeEchemFrame(EchemFrame, ResizeFrame)
            
        # Initialize plotter
        self.EchemFig = EchemFig(self.fig, self)
        self.MakePlotParamsFrame(PlotParamsFrame)
        self.MakeUpdateFrame(UpdateButtonFrame)
        self.MakePlotTypeFrame(PlotTypeFrame)
        
        self.fig2selection.trace('w', self.fig_opt_changed)
        self.fig2ptselection.trace('w', self.fig_opt_changed)
        
        # Norm traces
        self.CurrentUnits_.trace('w', self.fig_opt_changed)
        self.reference_var.trace('w', self.fig_opt_changed)
        self.PotentialUnits_.trace('w', self.fig_opt_changed)
        self.TimeUnits_.trace('w', self.fig_opt_changed)
        self.Colorbar_.trace('w', self.fig_opt_changed)
        self.Location_for_cbar.trace('w', self.fig_opt_changed)
        self.CurrentDensity_.trace('w', self.fig_opt_changed)
        self.curr_den_units.trace('w', self.fig_opt_changed)
        
        # Both
        self.Overlay_.trace('w', self.fig_opt_changed)
        self.plot_cmap.trace('w', self.fig_opt_changed)
        
        # EIS traces
        self.EIS_view_selection.trace('w', self.fig_opt_changed)
        self.ImpedanceUnits_.trace('w', self.fig_opt_changed)
        self.Scale_EIS.trace('w', self.fig_opt_changed)
        
        # Customize traces
        self.Legend_.trace('w', self.fig_opt_changed)
        self.location_for_legend.trace('w', self.fig_opt_changed)
        self.line_style.trace('w', self.fig_opt_changed)
        self.marker_style.trace('w', self.fig_opt_changed)
        self.Inset_.trace('w', self.fig_opt_changed)
        self.location_for_Inset.trace('w', self.fig_opt_changed)
        self.Ticks_Inset.trace('w', self.fig_opt_changed)
        self.spine_line_style_Inset.trace('w', self.fig_opt_changed)
        self.zoom_line_style_Inset.trace('w', self.fig_opt_changed)
        self.upperleft_line_style_Inset.trace('w', self.fig_opt_changed)
        self.upperright_line_style_Inset.trace('w', self.fig_opt_changed)
        self.lowerleft_line_style_Inset.trace('w', self.fig_opt_changed)
        self.lowerright_line_style_Inset.trace('w', self.fig_opt_changed)
        
        # Data control traces
        self.apply_notch_filter.trace('w', self.fig_opt_changed)
        self.AnalyzePeak_.trace('w', self.fig_opt_changed)
        self.baseline_analysis_line_style.trace('w', self.fig_opt_changed)
        
        # Partical import data traces
        self.start_after_option.trace('w', self.fig_opt_changed)
        self.end_before_option.trace('w', self.fig_opt_changed)
        
        # Send print messages to the console
        sys.stdout = PrintLogger(console)
        
    #################### END __init__ ####################
    
    ########## GUI CALLBACKS ##########   
    
    def openFolder(self):
        folder = filedialog.askdirectory()
        if folder == '':
            print('Error: Folder not chosen')
            return
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
    
    def macro_RS(self):
        if not hasattr(self, 'Macro_RS'):
            self.Popup_Generator = Popup_Generator()
        self.Popup_Generator.get(self, 'Macro_RS', None)
    
    def macro_RS_RDE(self):
        if not hasattr(self, 'Macro_RS_RDE'):
            self.Popup_Generator = Popup_Generator()
        self.Popup_Generator.get(self, 'Macro_RS_RDE', None)
    
    def micro_RS(self):
        if not hasattr(self, 'Micro_RS'):
            self.Popup_Generator = Popup_Generator()
        self.Popup_Generator.get(self, 'Micro_RS', None)
        
    def FT_data(self):
        # Get data from the lines
        x_data, y_data = [], []

        for line in self.EchemFig.ax.get_lines():
            x = line.get_xdata()
            y = line.get_ydata()
            if len(x) > 0 and len(y) > 0:
                x_data.append(x)
                y_data.append(y)

        if not x_data:
            print('Error: Cannot Fourier transform data. Load file first.')
            return
        
        # Keep as separate lists of arrays
        data = {"x": x_data, "y": y_data}
        
        if not hasattr(self, 'Popup_Generator'):
            self.Popup_Generator = Popup_Generator()
        self.Popup_Generator.get(self, 'FT_data', data)
    
    def copy_figure_to_clipboard(self):
        """Copy current matplotlib figure to system clipboard (cross-platform)."""
        buf = io.BytesIO()
        self.fig.savefig(buf, format="png", dpi=300, bbox_inches="tight")
        buf.seek(0)
        img_bytes = buf.getvalue()
        
        if sys.platform.startswith("win"):
            # --- Windows ---
            image = Image.open(io.BytesIO(img_bytes))
            output = io.BytesIO()
            image.convert("RGB").save(output, "BMP")
            data = output.getvalue()[14:]  # strip BMP header
            output.close()
    
            win32clipboard.OpenClipboard()
            win32clipboard.EmptyClipboard()
            win32clipboard.SetClipboardData(win32clipboard.CF_DIB, data)
            win32clipboard.CloseClipboard()
            print("✔️ Figure copied to clipboard 🪟")
        
        elif sys.platform == "darwin":
            # --- macOS ---
            p = subprocess.Popen(["pbcopy"], stdin=subprocess.PIPE)
            p.communicate(img_bytes)
            print("✔️ Figure copied to clipboard 🍎")
        
        elif sys.platform.startswith("linux"):
            # --- Linux ---
            try:
                p = subprocess.Popen(
                    ["xclip", "-selection", "clipboard", "-t", "image/png"],
                    stdin=subprocess.PIPE
                )
                p.communicate(img_bytes)
                print("✔️ Figure copied to clipboard (Linux)")
            except FileNotFoundError:
                print("⚠️ xclip not found. Install with: sudo apt install xclip")
        else:
            print("⚠️ Clipboard copy not supported on this platform.")
                
if __name__ == '__main__':
    root = Tk()
    gui = GUI()
    root.protocol("WM_DELETE_WINDOW", gui.exit_)
    root.mainloop()
