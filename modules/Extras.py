from tkinter import *
from tkinter.ttk import *
import matplotlib
import matplotlib.pyplot as plt
from matplotlib.backends.backend_tkagg import FigureCanvasTkAgg
from modules.Plotter import EchemFig
import numpy as np

class Make_Popup_GUI_Figure(EchemFig):
    
    def __init__(self, GUI):
        self.GUI = GUI
        self.fig = plt.Figure(figsize=(4,4), dpi=100, constrained_layout=True)
        self.fig.add_subplot(111)
        self.ax = self.fig.gca()
        self.dpi = 300
        self.EchemFig = EchemFig(self.fig, self)
        self.make_popup()
        self.fill_leftframe()
        self.draw()
        
    def make_popup(self):
        self.popup = Toplevel()
        self.leftframe  = Frame(self.popup)
        self.rightframe = Frame(self.popup)
        self.leftframe.grid(row=0, column=0)
        self.rightframe.grid(row=0, column=1)
        FigureCanvasTkAgg(self.fig, master=self.rightframe
                          ).get_tk_widget().grid(row=0, column=0)
        
    def fill_leftframe(self):
        pass
    
    def draw(self):
        pass
        
    def save(self):
        path = filedialog.asksaveasfilename(defaultextension='.png')
        if hasattr(self, 'dpi_field'):
            self.dpi = int(self.dpi_field.get())
        self.fig.savefig(path, dpi=self.dpi)
    
    def set_dpi(self, dpi):
        self.dpi = dpi
    
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
    
    def fill_leftframe(self):
        # Put relevant buttons in self.leftframe
        frame = self.leftframe
        
        if not hasattr(self, 'dpi_field'):
            # Will already have these attributes if it's being reinitialized
            self.dpi_field = StringVar(value=f'{self.dpi}')
            self.xlabel = StringVar(value='')
            self.ylabel = StringVar(value='')
            self.xticks = StringVar(value='')
            self.yticks = StringVar(value='')
            self.n_xticks = StringVar(value='')
            self.n_yticks = StringVar(value='')
            self.div_const = StringVar(value='1')
            self.draw_extras = IntVar(value=0)
        
        Label(frame, text='Fourier Transform Data    ').grid(row=0, column=0, columnspan=2)
        Button(frame, text='Redraw', command=self.redraw).grid(row=0, column=2, sticky=(W,E))
        
        Label(frame, text='dpi: ').grid(row=1, column=0, sticky=(W,E))
        Entry(frame, textvariable=self.dpi_field, width=4).grid(row=1, column=1,
                                                                sticky=(W,E))
        Button(frame, text='Save', command=self.save).grid(row=1, column=2, sticky=(W,E))
        
        
        Label(frame, text='X label: ').grid(row=2, column=0, 
                                                    columnspan=2, sticky=(W,E))
        Entry(frame, textvariable=self.xlabel, width=20).grid(
            row=2, column=1, columnspan=2, sticky=(W,E))
        
        Label(frame, text='Y label: ').grid(row=3, column=0, 
                                                    columnspan=2, sticky=(W,E))
        Entry(frame, textvariable=self.ylabel, width=20).grid(
            row=3, column=1, columnspan=2, sticky=(W,E))
        
        Label(frame, text='X ticks: ').grid(row=4, column=0, 
                                                    columnspan=2, sticky=(W,E))
        Entry(frame, textvariable=self.xticks, width=20).grid(
            row=4, column=1, columnspan=2, sticky=(W,E))
        
        Label(frame, text='Y ticks: ').grid(row=5, column=0, 
                                                    columnspan=2, sticky=(W,E))
        Entry(frame, textvariable=self.yticks, width=20).grid(
            row=5, column=1, columnspan=2, sticky=(W,E))
        
        Label(frame, text='# X ticks: ').grid(row=6, column=0, 
                                                    columnspan=2, sticky=(W,E))
        Entry(frame, textvariable=self.n_xticks, width=20).grid(
            row=6, column=1, columnspan=2, sticky=(W,E))
        
        Label(frame, text='# Y ticks: ').grid(row=7, column=0, 
                                                    columnspan=2, sticky=(W,E))
        Entry(frame, textvariable=self.n_yticks, width=20).grid(
            row=7, column=1, columnspan=2, sticky=(W,E))
        
        Label(frame, text='Div. const.: ').grid(row=8, column=0, 
                                                    columnspan=2, sticky=(W,E))
        Entry(frame, textvariable=self.div_const, width=20).grid(
            row=8, column=1, columnspan=2, sticky=(W,E))
        
        Checkbutton(frame, text='Draw annotations: ', variable=self.draw_extras
                    ).grid(row=9, column=0,columnspan=2, sticky=(W,E))
        
        pass
    
    def redraw(self):
        self.ax.clear()
        self.draw()
    
    
    def draw(self):
        x_array, y_array = self.data['x'], self.data['y']
        
        for x, y in x_array, y_array:
            y /= float(self.div_const.get())
            self.ax.plot(x,y, color='k')
            if self.draw_extras.get():
                self.draw_extra_artists()
            self.set_xlabel()
            self.set_ylabel()
            self.set_xticks()
            self.set_yticks()
            
            self.fig.canvas.draw()
            plt.pause(0.001)
    
    def divide_by_const(self, constant):
        pass
    
    
    def set_xlabel(self):
        xlabel = self.xlabel.get()
        self.ax.set_xlabel(f'{xlabel}')
    
    def set_ylabel(self):
        ylabel = self.ylabel.get()
        self.ax.set_ylabel(f'{ylabel}')
        
    def set_xticks(self):
        xticks = self.xticks.get()
        if xticks == '':
            return self.set_n_xticks()
        xticks = xticks.split(',')
        xticks = [float(n) for n in xticks]
        self.ax.set_xticks(xticks)
        
    def set_yticks(self):
        yticks = self.yticks.get()
        if yticks == '':
            return self.set_n_yticks()
        yticks = yticks.split(',')
        yticks = [float(n) for n in yticks]
        self.ax.set_yticks(yticks)
        
    def set_n_xticks(self):
        try:
            n_xticks = int(self.n_xticks.get())
        except:
            return
        self.ax.xaxis.set_major_locator(plt.MaxNLocator(n_xticks))
        
    def set_n_yticks(self):
        try:
            n_yticks = int(self.n_yticks.get())
        except:
            return
        self.ax.yaxis.set_major_locator(plt.MaxNLocator(n_yticks))
        
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