import numpy as np
from io import StringIO
import scipy.io
import os
import pandas as pd

'''
If multiple files to plot, add them to data folder and set Muti_files to True
'''
folder = r'Z:\Projects\Miguel\Raw data\2025\test\HEKA test\20250408.mat'
Multi_files = False
Overlay = None
Integrate = False

class extract_data():
    def __init__(self):
        return
        
    def read_file(self, folder, MF, GUI):
        # Local reference to GUI object
        self.GUI = GUI
        
        Multi_files = MF
        if Multi_files == False:
            file = folder
            Ts = None
            freq = None
            if file.endswith('.txt'):
                get_col = list(pd.read_csv(file, sep='\t', nrows=1).columns)
                if len(get_col) == 1:
                    get_col = list(pd.read_csv(file, sep=';', nrows=1).columns)
                # print(get_col)
                
                if get_col[0] == 'Index':
                    freq, real_Z, imag_Z, abs_Z, phase = extract.autolab_data_PEIS(file, 'Raw')
                    print(f'\nPloting Autolab EIS .txt file: {file.rsplit("/", 1)[-1]}')
                
                if get_col[0] == "Z' (Ω)":
                    freq, real_Z, imag_Z, abs_Z, phase = extract.autolab_data_PEIS(file, 'Fit')
                    print(f'\nPloting Autolab EIS fit .txt file: {file.rsplit("/", 1)[-1]}')
                
                if get_col[0] == "Potential applied (V)":
                    Ts, Vs, Is, sweeps = extract.autolab_data_norm(file)
                    print(f'\nPloting .txt file: {file.rsplit("/", 1)[-1]}')
                    NUM_SWEEPS = int(sweeps[-1])
                    # NUM_SWEEPS = 1
                    if NUM_SWEEPS == 0: #If no cycles, NUM_SWEEPS need to be set to 1
                        NUM_SWEEPS = 1
                    print(f'Number of cycles: {NUM_SWEEPS}')
                    
                if len(get_col) == 3:
                    Ts, Vs, Is = extract.bio_data_OVP(file)
                    print(f'\nPloting OCV .txt file: {file.rsplit("/", 1)[-1]}')
                    NUM_SWEEPS = 1
                
                if get_col[0] == 'freq/Hz':
                    if len(get_col) == 6:
                        freq, real_Z, imag_Z, abs_Z, phase = extract.bio_data_PEIS(file)
                        NUM_SWEEPS = 1
                    if len(get_col) == 7:
                        freq, real_Z, imag_Z, abs_Z, phase, sweeps = extract.bio_data_PEIS_cycles(file)
                        NUM_SWEEPS = int(sweeps[-1])
                        if NUM_SWEEPS == 0: #If no cycles, NUM_SWEEPS need to be set to 1
                            NUM_SWEEPS = 1
                    print(f'\nPloting EIS .txt file: {file.rsplit("/", 1)[-1]}')
                
                if len(get_col) == 4:
                    Ts, Vs, Is = extract.bio_data_no_cycle(file)
                    print(f'\nPloting .txt file: {file.rsplit("/", 1)[-1]}')
                    NUM_SWEEPS = 1
                
                if len(get_col) == 5:
                    Ts, Vs, Is, sweeps = extract.bio_data_norm(file)
                    print(f'\nPloting .txt file: {file.rsplit("/", 1)[-1]}')
                    NUM_SWEEPS = int(sweeps[-1])
                    if NUM_SWEEPS == 0: #If no cycles, NUM_SWEEPS need to be set to 1
                        NUM_SWEEPS = 1
                    print(f'Number of cycles: {NUM_SWEEPS}')
            
            if file.endswith('.mat'):
                while True:
                    try:
                        input_ = self.GUI.console_input('Input single series (integer) to plot\n>>\n')
                        series_to_plot = list(map(float, input_.split()))
                        if len(series_to_plot) > 1:
                            print('Warning: multiple series detected... only plotting first.')
                            break
                        elif len(series_to_plot) == 1:
                            break
                    except ValueError:
                        print(f'Invalid input "{input_}": Enter single integer value as series')
                
                series_to_plot = series_to_plot[0]
                try:
                    Ts, Vs, Is, sweeps = extract.matlab_iv_data(file, series_to_plot, False)
                except TypeError:
                    return None
                print(f'\nPloting HEKA Series {int(series_to_plot)} from .mat file: {file.rsplit("/", 1)[-1]}')
                NUM_SWEEPS = sweeps
                if NUM_SWEEPS == 0: #If no cycles, NUM_SWEEPS need to be set to 1
                    NUM_SWEEPS = 1
                print(f'Number of cycles: {NUM_SWEEPS}')
                
            if file.endswith('.asc'):
                Ts, Vs, Is, sweeps = extract.heka_data(file)
                print(f'\nPloting HEKA .asc file: {file.rsplit("/", 1)[-1]}')
                NUM_SWEEPS = sweeps[-1]
                if NUM_SWEEPS == 0: #If no cycles, NUM_SWEEPS need to be set to 1
                    NUM_SWEEPS = 1
                print(f'Number of cycles: {NUM_SWEEPS}')
                
            if file.endswith('.csv'):
                Ts, Vs, Is = extract.seccm_data(file)
                print(f'\nPloting .csv file: {file.rsplit("/", 1)[-1]}')
                NUM_SWEEPS = 1
                
            if Ts is not None:
                ##Split file by number of sweeps
                T, V, I = (np.array_split(Ts, NUM_SWEEPS),
                           np.array_split(Vs, NUM_SWEEPS),
                           np.array_split(Is, NUM_SWEEPS))
                extracted_data = [T, V, I, NUM_SWEEPS, 0, 'Norm']
                
            if freq is not None:
                freq, real_Z, imag_Z, abs_Z, phase = (np.array_split(freq, NUM_SWEEPS),
                                                      np.array_split(real_Z, NUM_SWEEPS),
                                                      np.array_split(imag_Z, NUM_SWEEPS),
                                                      np.array_split(abs_Z, NUM_SWEEPS),
                                                      np.array_split(phase, NUM_SWEEPS))
                extracted_data = [freq, real_Z, imag_Z, abs_Z, phase, NUM_SWEEPS, 0, 'EIS']
            
            return extracted_data, 0
        
        if Multi_files == True:
            extracted_data = []
            file_num = 0
            
            '''Biologic and Autolab files'''
            bio_files = sorted([os.path.join(folder, file)
                 for file in os.listdir(folder)
                 if file.endswith('.txt')])
            
            bio_file_max = len(bio_files)
            for f in bio_files:
                Ts = None
                freq = None
                ## Column names for data frame to distinguish exported .txt file
                get_col = list(pd.read_csv(f, sep='\t', nrows=1).columns)
                if len(get_col) == 1:
                    get_col = list(pd.read_csv(f, sep=';', nrows=1).columns)
                # print(get_col)
                
                if get_col[0] == 'Index':
                    freq, real_Z, imag_Z, abs_Z, phase = extract.autolab_data_PEIS(f, 'Raw')
                    print(f'\nPloting Autolab EIS .txt file: {f.rsplit("/", 1)[-1]}')
                
                if get_col[0] == "Z' (Ω)":
                    freq, real_Z, imag_Z, abs_Z, phase = extract.autolab_data_PEIS(f, 'Fit')
                    print(f'\nPloting Autolab EIS fit .txt file: {f.rsplit("/", 1)[-1]}')
                
                if get_col[0] == "Potential applied (V)":
                    Ts, Vs, Is, sweeps = extract.autolab_data_norm(f)
                    print(f'\nPloting .txt file {file_num+1}: {f.rsplit("/", 1)[-1]}')
                    NUM_SWEEPS = int(sweeps[-1])
                    # NUM_SWEEPS = 1
                    if NUM_SWEEPS == 0: #If no cycles, NUM_SWEEPS need to be set to 1
                        NUM_SWEEPS = 1
                    print(f'Number of cycles: {NUM_SWEEPS}')
                    
                if len(get_col) == 3:
                    Ts, Vs, Is = extract.bio_data_OVP(f)
                    NUM_SWEEPS = 1
                    print(f'\nPloting OCV .txt file {file_num+1}: {f.rsplit("/", 1)[-1]}')
                            
                if get_col[0] == 'freq/Hz':
                    if len(get_col) == 6:
                        freq, real_Z, imag_Z, abs_Z, phase = extract.bio_data_PEIS(f)
                        NUM_SWEEPS = 1
                    if len(get_col) == 7:
                        freq, real_Z, imag_Z, abs_Z, phase, sweeps = extract.bio_data_PEIS_cycles(f)
                        NUM_SWEEPS = int(sweeps[-1])
                        if NUM_SWEEPS == 0: #If no cycles, NUM_SWEEPS need to be set to 1
                            NUM_SWEEPS = 1
                    print(f'\nPloting EIS .txt file {file_num+1}: {f.rsplit("/", 1)[-1]}')
                    
                if len(get_col) == 4:
                    Ts, Vs, Is = extract.bio_data_no_cycle(f)
                    print(f'\nPloting .txt file {file_num+1}: {f.rsplit("/", 1)[-1]}')
                    NUM_SWEEPS = 1
                    print(f'Number of cycles: {NUM_SWEEPS}')
                
                if len(get_col) == 5:
                    Ts, Vs, Is, sweeps = extract.bio_data_norm(f)
                    print(f'\nPloting .txt file {file_num+1}: {f.rsplit("/", 1)[-1]}')
                    NUM_SWEEPS = int(sweeps[-1])
                    # NUM_SWEEPS = 1
                    if NUM_SWEEPS == 0: #If no cycles, NUM_SWEEPS need to be set to 1
                        NUM_SWEEPS = 1
                    print(f'Number of cycles: {NUM_SWEEPS}')
                
                if Ts is not None:
                    ##Split file by number of sweeps
                    T, V, I = (np.array_split(Ts, NUM_SWEEPS),
                               np.array_split(Vs, NUM_SWEEPS),
                               np.array_split(Is, NUM_SWEEPS))
                    extracted_data.append([T, V, I, NUM_SWEEPS, file_num, 'Norm'])
                    
                if freq is not None:
                    freq, real_Z, imag_Z, abs_Z, phase = (np.array_split(freq, NUM_SWEEPS),
                                                          np.array_split(real_Z, NUM_SWEEPS),
                                                          np.array_split(imag_Z, NUM_SWEEPS),
                                                          np.array_split(abs_Z, NUM_SWEEPS),
                                                          np.array_split(phase, NUM_SWEEPS))
                    extracted_data.append([freq, real_Z, imag_Z, abs_Z, phase, NUM_SWEEPS, file_num, 'EIS'])
                
                file_num += 1
        
            '''SECCM or any .csv files'''
            seccm_files = sorted([os.path.join(folder, file)
                 for file in os.listdir(folder)
                 if file.endswith('.csv')])
            
            seccm_file_max = len(seccm_files)
            for f in seccm_files:
                Ts, Vs, Is = extract.seccm_data(f)
                print(f'\nPloting .csv file {file_num+1}: {f.rsplit("/", 1)[-1]}')
                NUM_SWEEPS = 1
                ##Split file by number of sweeps
                T, V, I = (np.array_split(Ts, NUM_SWEEPS),
                           np.array_split(Vs, NUM_SWEEPS),
                           np.array_split(Is, NUM_SWEEPS))
                extracted_data.append([T, V, I, NUM_SWEEPS, file_num, 'Norm'])
                file_num += 1
        
            '''HEKA files'''
            heka_files = sorted([os.path.join(folder, file)
                 for file in os.listdir(folder)
                 if file.endswith('.mat') or file.endswith('.asc')])
            
            heka_file_max = 0
            for f in heka_files:
                if f.endswith('.mat'):
                    while True:
                        try:
                            input_ = self.GUI.console_input('Input groups to plot with spaces in between (Example:1 2 3 4) or type "all" for all groups\n>>\n')
                            if input_ == 'all':
                                groups = extract.matlab_iv_data(f, 1, True)
                                print(f'Extracting all groups of length {groups}')
                                series_to_plot = list(range(1, groups))
                                # print(len(series_to_plot))
                                break
                                
                            series_to_plot = list(map(float, input_.split()))
                            if len(series_to_plot) >= 1:
                                break
                        except ValueError:
                            print(f'Invalid input "{input_}"')
                    
                    heka_file_max = len(series_to_plot)
                    for plot_ in series_to_plot: 
                        try:
                            Ts, Vs, Is, sweeps = extract.matlab_iv_data(f, int(plot_), False)
                        except TypeError:
                            return None
                        
                        print(f'\nPloting HEKA Series {int(plot_)} from .mat file: {f.rsplit("/", 1)[-1]}')
                        NUM_SWEEPS = sweeps
                        if NUM_SWEEPS == 0: #If no cycles, NUM_SWEEPS need to be set to 1
                            NUM_SWEEPS = 1
                        print(f'Number of cycles: {NUM_SWEEPS}')
                        ##Split file by number of sweeps
                        T, V, I = (np.array_split(Ts, NUM_SWEEPS),
                                   np.array_split(Vs, NUM_SWEEPS),
                                   np.array_split(Is, NUM_SWEEPS))
                        extracted_data.append([T, V, I, NUM_SWEEPS, file_num, 'Norm'])
                        
                        file_num += 1
                    
                if f.endswith('.asc'):
                    Ts, Vs, Is, sweeps = extract.heka_data(f)
                    print(f'\nPloting HEKA .asc file: {f.rsplit("/", 1)[-1]}')
                    NUM_SWEEPS = sweeps[-1]
                    if NUM_SWEEPS == 0: #If no cycles, NUM_SWEEPS need to be set to 1
                        NUM_SWEEPS = 1
                    print(f'Number of cycles: {NUM_SWEEPS}')
                    ##Split file by number of sweeps
                    T, V, I = (np.array_split(Ts, NUM_SWEEPS),
                               np.array_split(Vs, NUM_SWEEPS),
                               np.array_split(Is, NUM_SWEEPS))
                    extracted_data.append([T, V, I, NUM_SWEEPS, file_num, 'Norm'])
                    
                    heka_file_max += 1
                    file_num += 1
                    
            file_max = bio_file_max + seccm_file_max + heka_file_max
            
            return extracted_data, file_max
    
class extract:
    
    def bio_data_OVP(file):
        df = pd.read_csv(file, names=('t', 'v'), skiprows=1, sep='\t')
        
        t = np.array(df['t'])
        v = np.array(df['v'])
        i = [0 for _ in range(len(t))]
        
        return t, v, i
    
    def autolab_data_PEIS(file, type_):
        if type_ == 'Raw':
            df = pd.read_csv(file, names=('Index', 'freq', 'r(Z)', 'im(Z)', '|Z|', 'phase', 't'), skiprows=1, sep=';')
            
        if type_ == 'Fit':
            df = pd.read_csv(file, names=('r(Z)', 'im(Z)', 'Error r(Z)', 'Error im(Z)', '|Z|', 'phase', 'freq', 'convergence', 'number of iterations', 'Chi-squared'), skiprows=1, sep=';')
        
        freq = np.array(df['freq'])
        real_Z = np.array(df['r(Z)'])
        imag_Z = np.array(df['im(Z)'])
        abs_Z = np.array(df['|Z|'])
        phase = np.array(df['phase'])
        
        return freq, real_Z, imag_Z, abs_Z, phase
    
    def bio_data_PEIS(file):
        df = pd.read_csv(file, names=('freq', 'r(Z)', 'im(Z)', '|Z|', 'phase'), skiprows=1, sep='\t')
        
        freq = np.array(df['freq'])
        real_Z = np.array(df['r(Z)'])
        imag_Z = np.array(df['im(Z)'])
        abs_Z = np.array(df['|Z|'])
        phase = np.array(df['phase'])
        
        return freq, real_Z, imag_Z, abs_Z, phase
    
    def bio_data_PEIS_cycles(file):
        df = pd.read_csv(file, names=('freq', 'r(Z)', 'im(Z)', '|Z|', 'phase', 'nc'), skiprows=1, sep='\t')
        
        freq = np.array(df['freq'])
        real_Z = np.array(df['r(Z)'])
        imag_Z = np.array(df['im(Z)'])
        abs_Z = np.array(df['|Z|'])
        phase = np.array(df['phase'])
        sweeps = np.array(df['nc'])
        
        return freq, real_Z, imag_Z, abs_Z, phase, sweeps
    
    def autolab_data_norm(file):
        df = pd.read_csv(file, names=('Potential applied', 't', 'i', 'v', 'cycle number', 'Index', 'Q+', 'Q-', 'current range'), skiprows=1, sep=';')
        
        t = np.array(df['t'])
        v = np.array(df['v'])
        i = np.array(df['i'])
        sweeps = np.array(df['cycle number'])
        
        return t, v, i, sweeps
    
    def bio_data_no_cycle(file):
        df = pd.read_csv(file, names=('t', 'v', 'i'), skiprows=1, sep='\t')
        
        t = np.array(df['t'])
        v = np.array(df['v'])
        i = np.array(df['i'])*1e-3            # Conversion to amps. i.e. data in mA
        
        return t, v, i
    
    def bio_data_norm(file):
        df = pd.read_csv(file, names=('t', 'v', 'i', 'cycle number'), skiprows=1, sep='\t')
        
        t = np.array(df['t'])
        v = np.array(df['v'])
        i = np.array(df['i'])*1e-3            # Conversion to amps. i.e. data in mA
        sweeps = np.array(df['cycle number'])
        
        return t, v, i, sweeps
    
    def seccm_data(file):
        df = pd.read_csv(file, names=('t', 'v', 'i'), skiprows=1,)
        
        t = np.array(df['t'])
        v = np.array(df['v'])
        i = np.array(df['i'])
        
        return t, v, i
    
    def heka_data(file):
        '''
        Parse PATCHMASTER-output csv files.
        
        Use StringIO to parse through file (fastest method I've found)
        Convert only floats to np arrays
        '''
        
        def isFloat(x):
            try: 
                float(x)
                return True
            except: 
                return False
        
        sweeps = []
        s = StringIO()
        with open(file, 'r') as f:
            for line in f:
                if isFloat(line.split(',')[0]):
                    # Check for index number
                    s.write(line)
                if line.startswith('Sweep'):
                    # print(line)
                    # print(line.split(','))
                    sweep_text = line.split(',')
                    sweep_text = sweep_text[0].split('_')
                    # print(sweep_text)
                    sweeps.append(int(sweep_text[-1]))
        
        if s.getvalue() != '':
            s.seek(0)
            array = np.genfromtxt(s, delimiter=',')
            array = array.T
            _, t, i, _, v = array
        
        return t, v, i, sweeps
    
    def matlab_iv_data(file, plot, extract_groups):
        '''
        Extracts I-V type data from a PATHCMASTER-generated .mat file.
        Assumes Trace 1 in PATCHMASTER is I (Current)
        Assumes Trace 2 in PATCHMASTER is V (Voltage)
        
        Returns times, voltages, currents
        
        Exporting binary .mat files from PATCHMASTER is much, much faster
        than exporting as csv files    
        '''
        d = scipy.io.loadmat(file)
        # print(d)
        groups = []
        series = []
        sweeps = []
        traces = []
        
        for key in d.keys():
            if not key.startswith('Trace'):
                continue
            root , group, season, sweep, trace = key.split('_')
            groups.append(int(group))
            series.append(int(season))
            sweeps.append(int(sweep))
            traces.append(int(trace))
        
        groups = np.array(groups)
        series = np.array(series)
        sweeps = np.array(sweeps)
        traces = np.array(traces)
        
        groups = groups[np.where(traces == 1)]
        series = series[np.where(traces == 1)]
        sweeps = sweeps[np.where(traces == 1)]
        traces = traces[np.where(traces == 1)]
        
        if extract_groups == True:
            return len(groups)
        
        if plot > len(groups):
            # print(len(groups))
            print(f'Warning: Cycle input out of bounds for data with length {len(groups)}')
            return None
        
        data = dict()
        count = 0
        # print(groups, series, sweeps, traces)
        while count < len(groups):
            
            T = np.array([])
            V = np.array([])
            I = np.array([])
            
            trace = 1
            key = f'{root}_{groups[count]}_{series[count]}_{sweeps[count]}_{trace}'
            arr = d[key]           # [ [t1, v1], [t2, v2], ...]
            arr = arr.transpose()  # [ [t1, t2, ...], [v1, v2, ...] ]
            # print(arr)
            ts, vals = arr
            
            if trace == traces[0]:
                T = np.append(T, ts)
            if trace == 1:
                I = np.append(I, vals)
            
            trace = 2
            key = f'{root}_{groups[count]}_{series[count]}_{sweeps[count]}_{trace}'
            arr = d[key]           # [ [t1, v1], [t2, v2], ...]
            arr = arr.transpose()  # [ [t1, t2, ...], [v1, v2, ...] ]
            # print(arr)
            ts, vals = arr
            
            if trace == 2:
                V = np.append(V, vals)
            
            if series[count] in data:
                data[series[count]].append([T, V, I]) # [ [[T1],[V1],[I1]], [[T2],[V2],[I2]], ...]
            else:
                data[series[count]] = list([[T, V, I]])
            count += 1
        
        t = []
        v = []
        i = []
        
        # print(len(data[plot]))
        # print(len(data[plot][0]))
        for cycles in data[plot]:
            t.append(cycles[0]) # [[T1], [T2], [T3], ...]
            v.append(cycles[1]) # [[V1], [V2], [V3], ...]
            i.append(cycles[2]) # [[I1], [I2], [I3], ...]
        # print(len(t))
        
        '''From: https://stackoverflow.com/questions/952914/how-do-i-make-a-flat-list-out-of-a-list-of-lists
        Flaten list of lists if multiple cycles'''
        
        if len(t)>1:
            t = [[tss for ts in t for tss in ts]]   # [[T1, T2, T3, ...]]
            v = [[vss for vs in v for vss in vs]]   # [[V1, V2, V3, ...]]
            i = [[iss for i_s in i for iss in i_s]] # [[I1, I2, I3, ...]]
        # print(len(t))
        
        t_return = t[0]
        v_return = v[0]
        i_return = i[0]
        
        # Option to import partical data set
        START_AFTER = False
        if type(START_AFTER) == int:    
            variable = t_return    # Choose variable to cut
            t_return = t_return[variable >= START_AFTER]
            v_return = v_return[variable >= START_AFTER]
            i_return = i_return[variable >= START_AFTER]
            print('Partial data extracted (Start {variable} at {START_AFTER})')
        
        END_BEFORE = False
        if type(END_BEFORE) == int:
            variable = t_return    # Choose variable to cut
            t_return = t_return[variable <= END_BEFORE]
            v_return = v_return[variable <= END_BEFORE]
            i_return = i_return[variable <= END_BEFORE]
            print('Partial data extracted (End {variable} at {END_BEFORE})')
        
        return t_return, v_return, i_return, len(data[plot])

    
if __name__ == '__main__':
    file = extract_data.read_file(folder, Multi_files)                          