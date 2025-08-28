import matplotlib as mpl
from matplotlib import pyplot as plt, cm, ticker as tk
from matplotlib.colors import LinearSegmentedColormap
from matplotlib.backends.backend_tkagg import FigureCanvasTkAgg
from mpl_toolkits.axes_grid1.anchored_artists import AnchoredSizeBar
from mpl_toolkits.axes_grid1 import make_axes_locatable
from matplotlib.colors import ListedColormap
from matplotlib.patches import Rectangle
import numpy as np
import os
from tkinter import *
from tkinter.ttk import *
from tkinter import filedialog
import sys
# Get the directory of the current script
script_dir = os.path.dirname(os.path.abspath(__file__))

# Add a directory relative to the script's location
sys.path.append(os.path.join(script_dir, ))

# from utils import Logger, nearest, threads

bool_map = {"true": True, 'True': True,
            "false": False, 'False': False}

time_unit_conv = {'hrs': 3600,
             'min': 60,
             's': 1,
             'ms': 1E-3,
             '\u03BCs': 1E-6}

potential_unit_conv = {'kV': 1000,
                       'V': 1,
                       'mV': 1E-3,
                       '\u03BCV': 1E-6}

current_unit_conv = {'fA': 1E-15,
             'pA': 1E-12,
             'nA': 1E-9,
             '\u03BCA': 1E-6,
             'mA': 1E-3,
             'A': 1}

Impedance_units_conv = {'G\u03A9': 1E9,
                        'M\u03A9': 1E6,
                        'k\u03A9': 1E3,
                        '\u03A9': 1}

def print_point(d:float,dec=0):
    '''
    Returns value as string with SI unit prefix
    
    e.g. unit_label(1e-9) --> '1 n'
         unit_label(7.7e-10): --> '770 p'
    '''
    inc_prefixes = ['k', 'M', 'G', 'T', 'P', 'E']
    dec_prefixes = ['m', 'µ', 'n', 'p', 'f', 'a']

    if d == 0:
        return f'0.0'

    degree = int(np.floor(np.log10(np.fabs(d)) / 3))

    prefix = ''

    if abs(degree) > 1:
        sign = degree / np.fabs(degree)
        if sign == 1:
            if degree - 1 < len(inc_prefixes):
                prefix = inc_prefixes[degree - 1]
            else:
                prefix = inc_prefixes[-1]
                degree = len(inc_prefixes)

        elif sign == -1:
            if -degree - 1 < len(dec_prefixes):
                prefix = dec_prefixes[-degree - 1]
            else:
                prefix = dec_prefixes[-1]
                degree = -len(dec_prefixes)

        scaled = float(d * pow(1000, -degree))

        s = f"{scaled:0.{dec}f}".rjust(4, ' ') + f" {prefix}"

    else:
        s = f"{d:0.4f}".rjust(4, ' ')
    return s


class EchemFig():
    def __init__(self, fig, GUI):
        # Local reference to GUI object
        self.GUI = GUI
        
        # Set up this figure, axes, and line to plot to
        self.fig = fig
        self.ax  = fig.gca()
        
        if not hasattr(self, "divider"):
            self.divider = make_axes_locatable(self.ax)

        self.ln, = self.ax.plot([0,1],[0,1])
        # self.load_style_file("Z:\Projects\Miguel\Spyder\style.mplstyle")
        
        # Keep track of what is currently being plotted
        self.DataPoint = None
        self.artists = []
        
        self.initialize()
        
        # Save original limits for reset
        self.original_xlim = self.ax.get_xlim()
        self.original_ylim = self.ax.get_ylim()
        
        # Rectangle patch (selection box)
        self.rect = Rectangle((0, 0), 0, 0, fill=False, color="red",
                              linewidth=1.5, linestyle="--")
        self.ax.add_patch(self.rect)
        self.rect.set_visible(False)

        # Event connections
        self.press_event = None
        self.cid_press = self.fig.canvas.mpl_connect("button_press_event", self.on_press)
        self.cid_release = self.fig.canvas.mpl_connect("button_release_event", self.on_release)
        self.cid_motion = self.fig.canvas.mpl_connect("motion_notify_event", self.on_motion)

    def on_press(self, event):
        if event.inaxes != self.ax:  # ignore clicks outside plot
            return
        self.press_event = event
        # Create new rectangle on the fly
        self.rect = Rectangle((event.xdata, event.ydata), 0, 0,
                              fill=False, color="k",
                              linewidth=1, linestyle="--")
        self.ax.add_patch(self.rect)
        self.rect.set_width(0)
        self.rect.set_height(0)
        self.fig.canvas.draw_idle()

    def on_motion(self, event):
        if self.press_event is None or event.inaxes != self.ax:
            return
        # Update rectangle as you drag
        x0, y0 = self.press_event.xdata, self.press_event.ydata
        x1, y1 = event.xdata, event.ydata
        self.rect.set_xy((min(x0, x1), min(y0, y1)))
        self.rect.set_width(abs(x1 - x0))
        self.rect.set_height(abs(y1 - y0))
        self.fig.canvas.draw_idle()

    def on_release(self, event):
        if self.press_event is None or event.inaxes != self.ax:
            return

        # Coordinates
        x0, y0 = self.press_event.xdata, self.press_event.ydata
        x1, y1 = event.xdata, event.ydata

        dx, dy = abs(x1 - x0), abs(y1 - y0)

        if dx < 1e-6 and dy < 1e-6:
            # Treat as simple click → print coords
            x, y = event.xdata, event.ydata
            print(f'({print_point(x, dec=3)}, {print_point(y, dec=3)})')
        else:
            # Treat as zoom box
            self.ax.set_xlim(min(x0, x1), max(x0, x1))
            self.ax.set_ylim(min(y0, y1), max(y0, y1))

        # Reset state
        self.press_event = None
        self.rect.set_visible(False)
        self.fig.canvas.draw_idle()

    def reset_zoom(self):
        """Restore the original limits"""
        self.ax.set_xlim(self.original_xlim)
        self.ax.set_ylim(self.original_ylim)
        self.fig.canvas.draw_idle()
        
    def initialize(self):
        '''
        Clear the figure
        '''
        self.ax.cla()

        self.ln, = self.ax.plot([],[])
        self.clear_artists()
        self.fig.tight_layout()
        
        if not hasattr(self, "style_dir") or self.style_dir is None:
            self.style_dir = None
        else:
            plt.style.use(self.style_dir)
            
        try:
            if hasattr(self, "cbar") and self.cbar:
                self.cbar.remove()
                self.cbar = None
                print('Colorbar removed')
            if hasattr(self, "cax") and self.cax:
                self.cax.remove()
                self.divider.remove()
                self.cax = None
                print('Divider removed')
        except Exception:
            pass
            
        self.fig.canvas.draw()
        # self.draw_artists()
    
            
    def set_data_from_file(self, extracted_data, file_max):
        self.file_max = file_max
        self.extracted_data = extracted_data
        # print('Data in Plotter')
    
    def update_plot(self):
        '''
        Display a given DataPoint by changing the selected 
        viewing option (from I vs V to I vs t, for example)
        
        Returns
        -------
        None.
        '''
        
        if not hasattr(self, "extracted_data") or self.extracted_data is None:
            # print("Error: No data loaded yet.")
            return
        
        width = int(self.GUI.plot_width.get())
        height = int(self.GUI.plot_height.get())
        self.fig.set_size_inches(width, height)
        
        # self.fig.tight_layout()
            
        overlay_val = self.GUI.Overlay_.get()
        # Convert string to boolean using the dictionary
        Overlay = bool_map.get(overlay_val.strip().lower(), False)
        
        if Overlay == False:
            self.initialize()
            file_to_plot = self.GUI.fig2ptselection.get()-1
            
            if self.extracted_data[file_to_plot][-1] == 'Norm':
                self.T, self.V, self.I, self.NUM_SWEEPS, self.file_num, self.file_type = self.extracted_data[file_to_plot]
                self.set_params_IV()
                
            if self.extracted_data[file_to_plot][-1] == 'EIS':
                self.freq, self.real_Z, self.imag_Z, self.abs_Z, self.phase, self.file_num, self.file_type = self.extracted_data[file_to_plot]
                self.set_params_EIS()
            print(f'Plotting file {file_to_plot+1}')
            return
        
        types_of_plots = []
        for file in range(self.file_max):
            type_ = self.extracted_data[file][-1]
            if type_ not in types_of_plots:
                types_of_plots.append(type_)
                
        if len(types_of_plots) == 2:
            # print(types_of_plots)
            print('\nError: Cannot Overlay EIS and Voltamperometric plots. Separate them in different folders!\n')
            return
        
        for file in range(self.file_max):
            if file == 0:
                self.initialize()
                
            if self.extracted_data[file][-1] == 'Norm':
                self.T, self.V, self.I, self.NUM_SWEEPS, self.file_num, self.file_type = self.extracted_data[file]
                self.set_params_IV()
                
            if self.extracted_data[file][-1] == 'EIS':
                self.freq, self.real_Z, self.imag_Z, self.abs_Z, self.phase, self.file_num, self.file_type = self.extracted_data[file]
                self.set_params_EIS()
        
        print('Plot updated')
    
    def set_params_IV(self):
        '''
        Grab display options from the GUI and initiate redrawing the plot.

        Returns
        -------
        None.
        '''
        IV_selection = self.GUI.fig2selection.get()
        
        overlay_val = self.GUI.Overlay_.get()
        # Convert string to boolean using the dictionary
        Overlay = bool_map.get(overlay_val.strip().lower(), False)
        
        # Unit and axis labels
        time_unit_val = self.GUI.TimeUnits_.get()
        ref_val = self.GUI.Reference_update.get()
        ref_unit_val = self.GUI.PotentialUnits_.get()
        current_unit_val = self.GUI.CurrentUnits_.get()
        IV_unit_params = [time_unit_val, ref_val, ref_unit_val, current_unit_val]
        
        # print('Params Set')
        self.plot_IV(IV_selection, IV_unit_params, Overlay)
        
             
    def plot_IV(self, IV_selection, IV_unit_params, Overlay):
        '''
        Function for drawing voltammetric (I-V, I-t, or V-t) data

        Parameters
        ----------
        IV_selection : string, display setting pulled from GUI

        Returns
        -------
        None.
        '''
        orig_pos = self.ax.get_position()
        
        if not hasattr(self, "T") or self.T is None:
            print("Error: No data loaded yet.")
            return
        
        d = {'t':self.T, 'V':self.V, 'I':self.I}
        
        time_units, ref_label, ref_units, current_units = IV_unit_params
        ylabel, xlabel = IV_selection.split(' vs ')
        
        count = 0
    
        # Decide how many colors are needed
        if self.file_max > 1 and Overlay == False and self.NUM_SWEEPS > 1:
            n_colors = self.NUM_SWEEPS
        else:
            n_colors = self.file_max if self.file_max > 1 else self.NUM_SWEEPS
        
        sm, colors, cmap = self.set_colormap_for_plot(n_colors)
        
        legend_labels = self.set_legend_labels(n_colors)
        
        marker_styles, marker_sizes, line_styles, line_sizes = self.set_line_markers_params(n_colors)
        
        x_shifts, y_shifts, cycles_to_plot = self.set_IV_cycles_to_plot_and_shifts(n_colors)
        
        while count < self.NUM_SWEEPS:
            yvals, xvals = d[ylabel][count], d[xlabel][count]
                
            if xlabel == 't':
                    xunits = time_unit_conv[time_units]
                    xvals = xvals/xunits
            
            if xlabel == 'V':
                    xunits = potential_unit_conv[ref_units]
                    xvals = xvals/xunits
            
            if ylabel == 'V':
                    yunits = potential_unit_conv[ref_units]
                    yvals = yvals/yunits
            
            if ylabel == 'I':
                    yunits = current_unit_conv[current_units]
                    yvals = yvals/yunits
            
            # Option to import partical data set
            start_mask = None
            end_mask = None
            skip_sweep = False
            
            Start_After = bool_map.get(self.GUI.start_after_option.get().strip().lower(), False)
            cutoff = self.GUI.start_after_float.get()
            if Start_After:
                if cutoff == '' and count == 0:
                    print('Error: input start after number')
                else:
                    key = self.GUI.start_after_var.get()
                    try:
                        ref_array = d[key][count]
                        if key == 't':
                            units = time_unit_conv[time_units]
                        if key == 'V':
                            units = potential_unit_conv[ref_units]
                        if key == 'I':
                            units = current_unit_conv[current_units]
                            
                        start_mask = ref_array/units > float(cutoff)
                    except ValueError:
                        pass
            
            End_Before = bool_map.get(self.GUI.end_before_option.get().strip().lower(), False)
            cutoff = self.GUI.end_before_float.get()
            if End_Before:
                if cutoff == '' and count == 0:
                    print('Error: input end before number')
                else:
                    key = self.GUI.end_before_var.get()
                    try:
                        ref_array = d[key][count]
                        if key == 't':
                            units = time_unit_conv[time_units]
                        if key == 'V':
                            units = potential_unit_conv[ref_units]
                        if key == 'I':
                            units = current_unit_conv[current_units]
                            
                        end_mask = ref_array/units < float(cutoff)
                    except ValueError:
                        pass
                    
            # ----- Apply masks -----
            if start_mask is not None and end_mask is not None:
                mask = start_mask & end_mask
            elif start_mask is not None:
                mask = start_mask
            elif end_mask is not None:
                mask = end_mask
            else:
                mask = None
            
            # ----- Apply mask safely -----
            if mask is not None:
                if len(mask) != len(xvals):  
                    print(f"Warning: skipping filtering sweep {count}, mask length mismatch.")
                else:
                    xvals = xvals[mask]
                    yvals = yvals[mask]
            
                    if len(xvals) == 0 or len(yvals) == 0:
                        print(f"Warning: no data left to plot after filtering sweep {count}. Skipping.")
                        skip_sweep = True
                    else:
                        skip_sweep = False
            
            if Overlay:
                color_idx = self.file_num
            else:
                color_idx = count
                
            count += 1
            
            # ✅ final check outside try/except to control plotting
            if skip_sweep:
                continue
            
            if color_idx in cycles_to_plot:
                self.ln, = self.ax.plot(xvals + x_shifts[color_idx],
                             yvals[:len(xvals)] + y_shifts[color_idx],
                             color=colors[color_idx],
                             marker = marker_styles[color_idx],
                             markersize = marker_sizes[color_idx],
                             linestyle = line_styles[color_idx],
                             linewidth = line_sizes[color_idx], 
                             label=legend_labels[color_idx])

            
            
        axis_labels = {'t': f'Time ({time_units})', 'V': f'E vs. {ref_label} ({ref_units})', 'I': f'Current ({current_units})'}
        # self.ax.set_box_aspect(1)
        self.ax.set_xlabel(axis_labels[xlabel])
        self.ax.set_ylabel(axis_labels[ylabel])
        self.ax.set_xscale('linear')
        self.update_plotlim()
        
        title = self.GUI.title_name.get()
        if title != '':
            try:
                self.ax.set_title(title)
            except ValueError:
                pass
        
        add_legend = bool_map.get(self.GUI.Legend_.get().strip().lower(), False)
        if add_legend == True:
            self.add_legend_to_plot()
                
        color_bar = bool_map.get(self.GUI.Colorbar_.get().strip().lower(), False)
        if color_bar == False:
            # Remove colorbar if it exists
            if hasattr(self, "cbar") and self.cbar:
                try:
                    self.cbar.remove()
                except:
                    pass
                self.cbar = None
                delattr(self, "cbar")
            if hasattr(self, "cax") and self.cax:
                try:
                    self.cax.remove()
                except:
                    pass
                self.cax = None
                delattr(self, "cax")
            
        if color_bar == True:
            # --- Create divider & cax only once ---
            location = self.GUI.Location_for_cbar.get()
            if location == 'top':
                orientation = 'horizontal'
            else:
                orientation = 'vertical'
                
            if not hasattr(self, "divider"):
                self.divider = make_axes_locatable(self.ax)
            if not hasattr(self, "cax"):
                self.cax = self.divider.append_axes(location, size="5%", pad=0.05)
                dummy_norm = mpl.colors.Normalize(vmin=1, vmax=1)
                dummy_sm = mpl.cm.ScalarMappable(norm=dummy_norm, cmap=cmap)
                self.cbar = self.fig.colorbar(dummy_sm, cax=self.cax, orientation=orientation)
            
            # location, label, aspect (ratio of long to short), fraction
            cbar_label = self.GUI.cbar_label.get()
            
            # --- Update colorbar instead of recreating ---
            self.cbar.update_normal(sm)
            self.cbar.set_label(cbar_label)
            
            # --- Get labels from GUI ---
            label_strs = self.GUI.labels_for_cbar.get().split(',')
            cbar_labels = [lbl.strip() for lbl in label_strs if lbl.strip() != '']
            try:
                num_ticks = int(self.GUI.ticknum_for_cbar.get().strip())
            except ValueError:
                num_ticks = None
                
            # Generate ticks
            if num_ticks and num_ticks > 1:
                cbar_ticks = np.linspace(1, n_colors, num_ticks)
            else:
                cbar_ticks = np.linspace(1, n_colors, n_colors)
            self.cbar.set_ticks(cbar_ticks)
            
            # If labels provided, must match number of ticks
            if cbar_labels and len(cbar_labels) != len(cbar_ticks):
                print("Warning: Number of labels does not match number of ticks.\nLabels are separated by commas\nIgnoring labels.")
                cbar_labels = []
            
            if cbar_labels:
                self.cbar.ax.set_yticklabels(cbar_labels)
            self.cax.set_position(self.ax.get_position())
            print(orientation)
            if orientation == 'horizontal':
                # Set the ticks position to the top
                self.cbar.ax.xaxis.set_ticks_position('top')
        
        self.ax.set_position(orig_pos)    
        self.draw_artists()
    
    def set_colormap_for_plot(self, n_colors):
        # Colormap option from GUI
        cmin, cmax = self.update_colormap()
        cmap = self.GUI.plot_cmap.get()
        color_params = [cmap, cmin, cmax]
        
        # Normalization maps sweep/file index -> colormap fraction
        norm = mpl.colors.Normalize(vmin=1, vmax=n_colors)
        
        # Colormap from your params
        if color_params[0] == 'Default':
            prop_cycle = mpl.rcParams['axes.prop_cycle']
            default_colors = prop_cycle.by_key()['color']
            # base_cmap = LinearSegmentedColormap.from_list("style_cycle_interp", default_colors, N=256)
            base_cmap = ListedColormap(default_colors, name="style_cycle_cmap")
        else:
            base_cmap = cm.get_cmap(color_params[0])
        
        if n_colors == 1:
            cmap = mpl.colors.ListedColormap([base_cmap(0)])
        else:
            cmap = self.truncate_colormap(base_cmap, color_params[1],
                                          color_params[2], n_colors)
        
        # ScalarMappable for both plotting & colorbar
        sm = mpl.cm.ScalarMappable(cmap=cmap, norm=norm)
        sm.set_array([])  # Matplotlib quirk
        colors = [sm.to_rgba(i) for i in range(1, n_colors+1)]
        
        return sm, colors, cmap
    
    def set_IV_cycles_to_plot_and_shifts(self, n_colors):
        # --- Get cycles to plot from GUI ---
        cycles_strs = self.GUI.cycles_to_plot.get().split(',')
        try:
            cycles_to_plot = [int(lbl.strip())-1 for lbl in cycles_strs if lbl.strip() != '']
        except Exception as e:
            print(f'Error: {e}\nCycles are integers separated by commas\nIgnoring cycles...')
        
        if len(cycles_to_plot)==0:
            cycles_to_plot = range(0, n_colors)
            
        # If cycles_to_plot provided, must match colors
        for num in cycles_to_plot:
            if num in range(0, n_colors):
                continue
            else:
                print("Warning: Cycles are out of bounds.\nCycles are integers separated by commas\nIgnoring cycles...")
                cycles_to_plot = range(0, n_colors)
                break
            
        # --- Get x shifts from GUI ---
        x_shifts_strs = self.GUI.x_axis_shifts.get().split(',')
        x_shifts = [float(lbl.strip()) for lbl in x_shifts_strs if lbl.strip() != '']
        
        if len(x_shifts)==0:
            x_shifts = [0 for _ in range(n_colors)]
            
        # If labels provided, must match number of colors
        if len(x_shifts)>0 and len(x_shifts) != n_colors:
            print("Warning: Number of x shifts does not match number of colors.\nShifts are separated by commas\nIgnoring shifts...")
            x_shifts = [0 for _ in range(n_colors)]
            
        # --- Get y shifts from GUI ---
        y_shifts_strs = self.GUI.y_axis_shifts.get().split(',')
        y_shifts = [float(lbl.strip()) for lbl in y_shifts_strs if lbl.strip() != '']
        
        if len(y_shifts)==0:
            y_shifts = [0 for _ in range(n_colors)]
            
        # If labels provided, must match number of colors
        if len(y_shifts)>0 and len(y_shifts) != n_colors:
            print("Warning: Number of y shifts does not match number of colors.\nShifts are separated by commas\nIgnoring shifts...")
            y_shifts = [0 for _ in range(n_colors)]
        
        return x_shifts, y_shifts, cycles_to_plot
    
    def set_legend_labels(self, n_colors):
        # --- Get labels from GUI ---
        label_strs = self.GUI.labels_for_legend.get().split(',')
        legend_labels = [lbl.strip() for lbl in label_strs if lbl.strip() != '']
        
        if len(legend_labels)==0:
            legend_labels = range(1, n_colors+1)
            
        # If labels provided, must match number of colors
        if len(legend_labels)>0 and len(legend_labels) != n_colors:
            print("Warning: Number of labels does not match number of colors.\nLabels are separated by commas\nIgnoring labels...")
            legend_labels = range(1, n_colors+1)
        
        return legend_labels
            
    def truncate_colormap(self, cmap, minval, maxval, n):
        """Return a truncated copy of a colormap"""
        new_cmap = LinearSegmentedColormap.from_list(
            f"trunc({cmap.name},{minval:.2f},{maxval:.2f})",
            cmap(np.linspace(minval, maxval, n)))
        return new_cmap

    def update_plotlim(self):
        '''
        Set Axis limits and tick multiples for IV plots
        '''
        # Grab values from GUI ('' if empty)
        x_multiple = self.GUI.x_axis_tickmultiple.get()
        xmin = self.GUI.echem_xminval.get()
        xmax = self.GUI.echem_xmaxval.get()
        
        y_multiple = self.GUI.y_axis_tickmultiple.get()
        ymin = self.GUI.echem_yminval.get()
        ymax = self.GUI.echem_ymaxval.get()
        # lim = (xmin, xmax), (ymin, ymax)
        
        x_minor_multiple = self.GUI.x_axis_minor_tickmultiple.get()
        y_minor_multiple = self.GUI.y_axis_minor_tickmultiple.get()
        
        xlim = (xmin, xmax)
        ylim = (ymin, ymax)
        
        # Update x-axis limits if both are provided
        if all(x != '' for x in xlim):
            try:
                self.ax.set_xlim(float(xlim[0]), float(xlim[1]))
            except Exception as e:
                print(f'Skiping setting x axis limit because of error: {e}')
                pass  # Could log or print warning
        
        # Update y-axis limits if both are provided
        if all(y != '' for y in ylim):
            try:
                self.ax.set_ylim(float(ylim[0]), float(ylim[1]))
            except Exception as e:
                print(f'Skiping setting y axis limit because of error: {e}')
                pass
        
        # Update x-axis tick multiple if provided
        if x_multiple != '':
            try:
                step = float(x_multiple)
                xmin, xmax = self.ax.get_xlim()
                n_ticks = abs((xmax - xmin) / step)
                if n_ticks > 100:
                    print(f"Error ⚠️ Too many major x-ticks ({int(n_ticks)}). Using AutoLocator instead...")
                    self.ax.xaxis.set_major_locator(tk.AutoLocator())
                else:
                    self.ax.xaxis.set_major_locator(tk.MultipleLocator(float(x_multiple)))
            except Exception as e:
                print(f'Skiping setting x axis major tick multiple because of error: {e}')
                pass
        
        # Update y-axis tick multiple if provided
        if y_multiple != '':
            try:
                step = float(y_multiple)
                ymin, ymax = self.ax.get_ylim()
                n_ticks = abs((ymax - ymin) / step)
                if n_ticks > 100:
                    print(f"Error ⚠️ Too many major y-ticks ({int(n_ticks)}). Using AutoLocator instead...")
                    self.ax.yaxis.set_major_locator(tk.AutoLocator())
                else:
                    self.ax.yaxis.set_major_locator(tk.MultipleLocator(float(y_multiple)))
            except Exception as e:
                print(f'Skiping setting y axis major tick multiple because of error: {e}')
                pass
        
        # Update x-axis minor tick multiple if provided
        if x_minor_multiple != '':
            try:
                step = float(x_minor_multiple)
                xmin, xmax = self.ax.get_xlim()
                n_ticks = abs((xmax - xmin) / step)
                if n_ticks > 100:
                    print(f"Error ⚠️ Too many minor x-ticks ({int(n_ticks)}). Using AutoLocator instead...")
                    self.ax.xaxis.set_minor_locator(tk.AutoLocator())
                else:
                    self.ax.xaxis.set_minor_locator(tk.MultipleLocator(float(x_minor_multiple)))
            except Exception as e:
                print(f'Skiping setting x axis minor tick multiple because of error: {e}')
                pass
        
        # Update y-axis minor tick multiple if provided
        if y_minor_multiple != '':
            try:
                step = float(y_multiple)
                ymin, ymax = self.ax.get_ylim()
                n_ticks = abs((ymax - ymin) / step)
                if n_ticks > 100:
                    print(f"Error ⚠️ Too many minor y-ticks ({int(n_ticks)}). Using AutoLocator instead...")
                    self.ax.yaxis.set_minor_locator(tk.AutoLocator())
                else:
                    self.ax.yaxis.set_minor_locator(tk.MultipleLocator(float(y_minor_multiple)))
            except Exception as e:
                print(f'Skiping setting y axis minor tick multiple because of error: {e}')
                pass
        
        self.original_xlim = self.ax.get_xlim()
        self.original_ylim = self.ax.get_ylim()
    
    def add_legend_to_plot(self):
        location = self.GUI.location_for_legend.get()
        handle_length = float(self.GUI.handle_length_legend.get())
        fontsize_ = self.GUI.fontsize_for_legend.get()
        if fontsize_ == '':
            fontsize_ = None
            
        self.ax.legend(fontsize=fontsize_,
                  labelcolor= 'k',
                  handlelength = handle_length,
                  loc = location)
    
    def set_line_markers_params(self, n_colors):
        '''Get line and marker styles from GUI'''
        marker_style = self.GUI.marker_style.get()
        marker_size = self.GUI.marker_width.get()
        line_style = self.GUI.line_style.get()
        line_size = self.GUI.line_width.get()
        
        # --- Get manual styles from GUI ---
        manual_ls_strs = self.GUI.manual_line_styles.get().split(',')
        manual_ls = [lbl.strip() for lbl in manual_ls_strs if lbl.strip() != '']
        
        manual_ms_strs = self.GUI.manual_marker_styles.get().split(',')
        manual_ms = [lbl.strip() for lbl in manual_ms_strs if lbl.strip() != '']
        
        # If no line styles provided
        if len(manual_ls)==0:
            line_style = [line_style for _ in range(n_colors)]
        # If line styles provided, must match number of colors
        elif len(manual_ls)>0 and len(manual_ls) != n_colors:
            print("Warning: Number of line styles does not match number of colors.\nLine styles are separated by commas\nIgnoring Manual LS...")
            line_style = [line_style for _ in range(n_colors)]
        else:
            line_style = manual_ls
        
        # If no marker styles provided
        if len(manual_ms)==0:
            marker_style = [marker_style for _ in range(n_colors)]
        # If marker styles provided, must match number of colors
        elif len(manual_ms)>0 and len(manual_ms) != n_colors:
            print("Warning: Number of marker styles does not match number of colors.\nMarker styles are separated by commas\nIgnoring Manual MS...")
            marker_style = [marker_style for _ in range(n_colors)]
        else:
            marker_style = manual_ms
        
        # Marker and Line size are always automatic
        marker_size = [marker_size for _ in range(n_colors)]
        line_size = [line_size for _ in range(n_colors)]
        
        return marker_style, marker_size, line_style, line_size
    
    def set_params_EIS(self):
        '''
        Grab display options from the GUI and initiate redrawing the plot.

        Returns
        -------
        None.
        '''
        
        EIS_selection = self.GUI.EIS_view_selection.get()
        
        overlay_val = self.GUI.Overlay_.get()
        # Convert string to boolean using the dictionary
        Overlay = bool_map.get(overlay_val.strip().lower(), False)
        
        # Unit and axis labels
        impedance_unit_val = self.GUI.ImpedanceUnits_.get()
        
        # Colormap option
        cmin, cmax = self.update_colormap()
        cmap = self.GUI.plot_cmap.get()
        color_params = [cmap, cmin, cmax]
        
        # print('Params Set')
        if EIS_selection == 'Nyquist':
            self.plot_Nyquist(Overlay, color_params, impedance_unit_val)
        if EIS_selection == '|Z| Bode':
            self.plot_Bode('Z', Overlay, color_params, impedance_unit_val)
        if EIS_selection == 'Phase Bode':
            self.plot_Bode('Phase', Overlay, color_params, impedance_unit_val)
        
    def plot_Nyquist(self, Overlay, color_params, impedance_units):
        '''
        Make a Nyquist plot.
        '''
        
        count = 0
        self.NUM_SWEEPS = 1
        
        # Decide how many colors are needed
        if self.file_max > 1 and Overlay == False and self.NUM_SWEEPS > 1:
            n_colors = self.NUM_SWEEPS
        else:
            n_colors = self.file_max if self.file_max > 1 else self.NUM_SWEEPS
        
        sm, colors, cmap = self.set_colormap_for_plot(n_colors)
        
        legend_labels = self.set_legend_labels(n_colors)
        
        marker_styles, marker_sizes, line_styles, line_sizes = self.set_line_markers_params(n_colors)
        
        x_shifts, y_shifts = self.set_EIS_shifts(n_colors)
        
        units = Impedance_units_conv[impedance_units]
        
        if Overlay:
            color_idx = self.file_num
        else:
            color_idx = count
            
        while count < self.NUM_SWEEPS:
            self.ax.plot(self.real_Z/units  + x_shifts[color_idx],
                         self.imag_Z/units  + y_shifts[color_idx],
                         color=colors[color_idx],
                         marker = marker_styles[color_idx],
                         markersize = marker_sizes[color_idx],
                         linestyle = line_styles[color_idx],
                         linewidth = line_sizes[color_idx], 
                         label = legend_labels[color_idx])
            count += 1
        
        self.ax.set_xscale(self.GUI.Scale_EIS.get())
        self.ax.set_xlabel(f"Z ' ({impedance_units})")
        self.ax.set_ylabel(f"- Z '' ({impedance_units})")
        self.update_plotlim_EIS()
        
        title = self.GUI.title_name_EIS.get()
        if title != '':
            try:
                self.ax.set_title(title)
            except ValueError:
                pass
        
        add_legend = bool_map.get(self.GUI.Legend_.get().strip().lower(), False)
        if add_legend == True:
            self.add_legend_to_plot()
            
        self.draw_artists()
    
    def plot_Bode(self, option, Overlay, color_params, impedance_units):
        '''
        Make a Bode plot
        '''
        
        count = 0
        self.NUM_SWEEPS = 1
        
        # Decide how many colors are needed
        if self.file_max > 1 and Overlay == False and self.NUM_SWEEPS > 1:
            n_colors = self.NUM_SWEEPS
        else:
            n_colors = self.file_max if self.file_max > 1 else self.NUM_SWEEPS
        
        sm, colors, cmap = self.set_colormap_for_plot(n_colors)
        
        legend_labels = self.set_legend_labels(n_colors)
        
        marker_styles, marker_sizes, line_styles, line_sizes = self.set_line_markers_params(n_colors) 
        
        x_shifts, y_shifts = self.set_EIS_shifts(n_colors)

        units = Impedance_units_conv[impedance_units]
        
        if Overlay:
            color_idx = self.file_num
        else:
            color_idx = count
            
        while count < self.NUM_SWEEPS:
            if option == 'Z':
                self.ax.plot(self.freq  + x_shifts[color_idx],
                             np.log10(self.abs_Z)/units + y_shifts[color_idx],
                             color=colors[color_idx],
                             marker = marker_styles[color_idx],
                             markersize = marker_sizes[color_idx],
                             linestyle = line_styles[color_idx],
                             linewidth = line_sizes[color_idx], 
                             label=legend_labels[color_idx])
                ylabel = f'log(|Z|) ({impedance_units})'
                
            elif option == 'Phase':
                self.ax.plot(self.freq  + x_shifts[color_idx], 
                             self.phase + y_shifts[color_idx],
                             color=colors[color_idx],
                             marker = marker_styles[color_idx],
                             markersize = marker_sizes[color_idx],
                             linestyle = line_styles[color_idx],
                             linewidth = line_sizes[color_idx], 
                             label=legend_labels[color_idx])
                ylabel = r'Phase ($\degree$)'
            
            count += 1
        
        self.ax.set_xscale(self.GUI.Scale_EIS.get())
        self.ax.set_xlabel('Frequency (Hz)')
        self.ax.set_ylabel(ylabel)
        self.update_plotlim_EIS()
        
        title = self.GUI.title_name_EIS.get()
        if title != '':
            try:
                self.ax.set_title(title)
            except ValueError:
                pass
            
        add_legend = bool_map.get(self.GUI.Legend_.get().strip().lower(), False)
        if add_legend == True:
            self.add_legend_to_plot()
            
        self.draw_artists()
    
    def update_plotlim_EIS(self):
        '''
        Set Axis limits and tick multiples for EIS plots
        '''
        x_multiple = self.GUI.x_axis_tickmultiple_EIS.get()
        xmin = self.GUI.echem_xminval_EIS.get()
        xmax = self.GUI.echem_xmaxval_EIS.get()
        
        y_multiple = self.GUI.y_axis_tickmultiple_EIS.get()
        ymin = self.GUI.echem_yminval_EIS.get()
        ymax = self.GUI.echem_ymaxval_EIS.get()
        # lim = (xmin, xmax), (ymin, ymax)
        
        x_minor_multiple = self.GUI.x_axis_minor_tickmultiple_EIS.get()
        y_minor_multiple = self.GUI.y_axis_minor_tickmultiple_EIS.get()
        
        xlim = (xmin, xmax)
        ylim = (ymin, ymax)
        
        # Update x-axis limits if both are provided
        if all(x != '' for x in xlim):
            try:
                self.ax.set_xlim(float(xlim[0]), float(xlim[1]))
            except Exception as e:
                print(f'Skiping setting x axis limit because of error: {e}')
                pass  # Could log or print warning
        
        # Update y-axis limits if both are provided
        if all(y != '' for y in ylim):
            try:
                self.ax.set_ylim(float(ylim[0]), float(ylim[1]))
            except Exception as e:
                print(f'Skiping setting y axis limit because of error: {e}')
                pass
        
        # Update x-axis tick multiple if provided
        if x_multiple != '':
            try:
                step = float(x_multiple)
                xmin, xmax = self.ax.get_xlim()
                n_ticks = abs((xmax - xmin) / step)
                if n_ticks > 100:
                    print(f"Error ⚠️ Too many major x-ticks ({int(n_ticks)}). Using AutoLocator instead...")
                    self.ax.xaxis.set_major_locator(tk.AutoLocator())
                else:
                    self.ax.xaxis.set_major_locator(tk.MultipleLocator(float(x_multiple)))
            except Exception as e:
                print(f'Skiping setting x axis major tick multiple because of error: {e}')
                pass
        
        # Update y-axis tick multiple if provided
        if y_multiple != '':
            try:
                step = float(y_multiple)
                ymin, ymax = self.ax.get_ylim()
                n_ticks = abs((ymax - ymin) / step)
                if n_ticks > 100:
                    print(f"Error ⚠️ Too many major y-ticks ({int(n_ticks)}). Using AutoLocator instead...")
                    self.ax.yaxis.set_major_locator(tk.AutoLocator())
                else:
                    self.ax.yaxis.set_major_locator(tk.MultipleLocator(float(y_multiple)))
            except Exception as e:
                print(f'Skiping setting y axis major tick multiple because of error: {e}')
                pass
        
        # Update x-axis minor tick multiple if provided
        if x_minor_multiple != '':
            try:
                step = float(x_minor_multiple)
                xmin, xmax = self.ax.get_xlim()
                n_ticks = abs((xmax - xmin) / step)
                if n_ticks > 100:
                    print(f"Error ⚠️ Too many minor x-ticks ({int(n_ticks)}). Using AutoLocator instead...")
                    self.ax.xaxis.set_minor_locator(tk.AutoLocator())
                else:
                    self.ax.xaxis.set_minor_locator(tk.MultipleLocator(float(x_minor_multiple)))
            except Exception as e:
                print(f'Skiping setting x axis minor tick multiple because of error: {e}')
                pass
        
        # Update y-axis minor tick multiple if provided
        if y_minor_multiple != '':
            try:
                step = float(y_multiple)
                ymin, ymax = self.ax.get_ylim()
                n_ticks = abs((ymax - ymin) / step)
                if n_ticks > 100:
                    print(f"Error ⚠️ Too many minor y-ticks ({int(n_ticks)}). Using AutoLocator instead...")
                    self.ax.yaxis.set_minor_locator(tk.AutoLocator())
                else:
                    self.ax.yaxis.set_minor_locator(tk.MultipleLocator(float(y_minor_multiple)))
            except Exception as e:
                print(f'Skiping setting y axis minor tick multiple because of error: {e}')
                pass
        
        self.original_xlim = self.ax.get_xlim()
        self.original_ylim = self.ax.get_ylim()
    
    def set_EIS_shifts(self, n_colors):
        # --- Get x shifts from GUI ---
        x_shifts_strs = self.GUI.x_axis_shifts_EIS.get().split(',')
        x_shifts = [float(lbl.strip()) for lbl in x_shifts_strs if lbl.strip() != '']
        
        if len(x_shifts)==0:
            x_shifts = [0 for _ in range(n_colors)]
            
        # If labels provided, must match number of colors
        if len(x_shifts)>0 and len(x_shifts) != n_colors:
            print("Warning: Number of x shifts does not match number of colors.\nShifts are separated by commas\nIgnoring shifts...")
            x_shifts = [0 for _ in range(n_colors)]
            
        # --- Get y shifts from GUI ---
        y_shifts_strs = self.GUI.y_axis_shifts_EIS.get().split(',')
        y_shifts = [float(lbl.strip()) for lbl in y_shifts_strs if lbl.strip() != '']
        
        if len(y_shifts)==0:
            y_shifts = [0 for _ in range(n_colors)]
            
        # If labels provided, must match number of colors
        if len(y_shifts)>0 and len(y_shifts) != n_colors:
            print("Warning: Number of y shifts does not match number of colors.\nShifts are separated by commas\nIgnoring shifts...")
            y_shifts = [0 for _ in range(n_colors)]
        
        return x_shifts, y_shifts
    
    def clear_artists(self):
        '''
        Removes all artists (except for self.ln) from the figure.

        Returns
        -------
        None.
        '''
        if len(self.artists) != 0:
            for artist in self.artists:
                try:
                    artist.remove()
                except:
                    # Artist may have never been drawn
                    pass
        self.artists = []
        
    
    def update_colormap(self):
        # cmap = self.GUI.plot_cmap.get()
        
        cmin = float(self.GUI.echem_cmap_minval.get())
        cmax = float(self.GUI.echem_cmap_maxval.get())
        
        # print('Color map updated')
        
        if (cmin < 0) or (cmin > 1) or (cmax < 0) or (cmax > 1):
            print('\nInvalid input! Min and max must be between 0 and 1.\n')
            return 0, 1
        
        return cmin, cmax
    
    def load_style_file(self, style_dir):
        self.style_dir = style_dir
        print('Style loaded')
        print(self.style_dir)
        
    def draw_artists(self):
        '''
        Draw all artists in the figure and draw the canvas.

        Returns
        -------
        None.
        '''
        self.ax.draw_artist(self.ln)
        for artist in self.artists:
            try:
                self.ax.add_artist(artist)
                self.ax.draw_artist(artist)
                self.fig.canvas.draw_idle()
            except Exception as e:
                '''If reloading an old data file, it fails to draw
                previously-generated artists. Catch that and remove the
                artist. (Can be regenerated by running the same analysis
                function)'''
                # print(f'Error drawing artist: {artist=}')
                # print(e)
                artist.remove()
                self.artists.remove(artist)  

        self.fig.tight_layout()
        self.fig.canvas.draw_idle()    
    
    def save(self):
        export_format = self.Export_format.get()
        path = filedialog.asksaveasfilename(defaultextension=export_format)
        if path == '':
            print('Error: plot not saved')
            return
        
        transparency = bool_map.get(self.Transparency.get().strip().lower(), False)
        dpi_val = int(self.dpi.get())
        
        self.fig.savefig(path, transparent = transparency, dpi=dpi_val, format=export_format)
        print(f'Plot exported as {export_format}')
        
    
    
    


