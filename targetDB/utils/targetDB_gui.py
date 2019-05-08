#!/usr/bin/env python
import re, sys
import tkinter as tk
from tkinter.scrolledtext import ScrolledText
from tkinter.messagebox import showerror, askokcancel, showinfo
from tkinter.filedialog import askopenfilename, askdirectory
import matplotlib.pyplot as plt
from matplotlib.backends.backend_tkagg import FigureCanvasTkAgg
import numpy as np
import configparser
from pathlib import Path

from targetDB import druggability_report as dr
from targetDB import target_descriptors as td
from targetDB.utils import gene2id as g2id


class ScrolledFrame(tk.Frame):
    def __init__(self, parent, vertical=True, horizontal=False):
        super().__init__(parent)

        # canvas for inner frame
        self._canvas = tk.Canvas(self)
        self._canvas.grid(row=0, column=0, sticky='news')  # changed

        # create right scrollbar and connect to canvas Y
        self._vertical_bar = tk.Scrollbar(self, orient='vertical', command=self._canvas.yview)
        if vertical:
            self._vertical_bar.grid(row=0, column=1, sticky='ns')
        self._canvas.configure(yscrollcommand=self._vertical_bar.set)

        # create bottom scrollbar and connect to canvas X
        self._horizontal_bar = tk.Scrollbar(self, orient='horizontal', command=self._canvas.xview)
        if horizontal:
            self._horizontal_bar.grid(row=1, column=0, sticky='we')
        self._canvas.configure(xscrollcommand=self._horizontal_bar.set)

        # inner frame for widgets
        self.inner = tk.Frame(self._canvas, bg='red')
        self._window = self._canvas.create_window((0, 0), window=self.inner, anchor='nw')

        # autoresize inner frame
        self.columnconfigure(0, weight=1)  # changed
        self.rowconfigure(0, weight=1)  # changed

        # resize when configure changed
        self.inner.bind('<Configure>', self.resize)
        self._canvas.bind('<Configure>', self.frame_width)

    def frame_width(self, event):
        # resize inner frame to canvas size
        canvas_width = event.width
        self._canvas.itemconfig(self._window, width=canvas_width)

    def resize(self, event=None):
        self._canvas.configure(scrollregion=self._canvas.bbox('all'))


class get_mpo_coeff_gui:
    def __init__(self):
        self.coeff = None
        self.master = tk.Toplevel()
        self.create_window()

        while True:
            self.master.update()
            if self.over:
                self.master.destroy()
                break

        while sum(self.coeff.values()) == 0:
            showerror('Error', 'Please set at least one coefficient to a value different from 0')
            self.master = tk.Toplevel()
            self.create_window()

            while True:
                self.master.update()
                if self.over:
                    self.master.destroy()
                    break

    def create_window(self):
        headers = ['Structural\ninformation', 'Structural\nDruggability', 'Chemistry', 'Biology', 'Diseases\nLinks',
                   'Genetic\nLinks', 'Literature\nInformation', 'Safety']
        colors = ['#95d0fc', '#0485d1', '#b790d4', '#87ae73', '#fec615', '#fb7d07', '#95a3a6', '#ff474c']
        tk.Label(self.master, height=2, text='Please input the desired weigth for the MPO score component',
                 font='bold').grid(row=0, columnspan=8)
        for i in range(len(headers)):
            tk.Label(self.master, height=2, text=headers[i], font='bold', bg=colors[i]).grid(row=2, column=i,
                                                                                             sticky=tk.W + tk.E)
        tk.Button(self.master, text='Validate', command=self.get_values, height=2, width=60, fg='red', bg='white').grid(
            columnspan=8, row=1)
        self.structural_info = tk.Scale(self.master, from_=200, to=-200, orient=tk.VERTICAL, length=300,
                                        tickinterval=50,
                                        bg='#95d0fc', font='bold')
        self.structural_info.set(100)
        self.structural_info.grid(row=3, column=0, sticky=tk.W + tk.E)
        self.structural_drug = tk.Scale(self.master, from_=200, to=-200, orient=tk.VERTICAL, length=300,
                                        tickinterval=50,
                                        bg='#0485d1', font='bold', fg='white')
        self.structural_drug.set(100)
        self.structural_drug.grid(row=3, column=1, sticky=tk.W + tk.E)
        self.chemistry = tk.Scale(self.master, from_=200, to=-200, orient=tk.VERTICAL, length=300, tickinterval=50,
                                  bg='#b790d4',
                                  font='bold')
        self.chemistry.set(100)
        self.chemistry.grid(row=3, column=2, sticky=tk.W + tk.E)
        self.biology = tk.Scale(self.master, from_=200, to=-200, orient=tk.VERTICAL, length=300, tickinterval=50,
                                bg='#87ae73',
                                font='bold')
        self.biology.set(100)
        self.biology.grid(row=3, column=3, sticky=tk.W + tk.E)
        self.disease = tk.Scale(self.master, from_=200, to=-200, orient=tk.VERTICAL, length=300, tickinterval=50,
                                bg='#fec615',
                                font='bold')
        self.disease.set(100)
        self.disease.grid(row=3, column=4, sticky=tk.W + tk.E)
        self.genetic = tk.Scale(self.master, from_=200, to=-200, orient=tk.VERTICAL, length=300, tickinterval=50,
                                bg='#fb7d07',
                                font='bold')
        self.genetic.set(100)
        self.genetic.grid(row=3, column=5, sticky=tk.W + tk.E)
        self.information = tk.Scale(self.master, from_=200, to=-200, orient=tk.VERTICAL, length=300, tickinterval=50,
                                    bg='#95a3a6', font='bold', fg='white')
        self.information.set(100)
        self.information.grid(row=3, column=6, sticky=tk.W + tk.E)
        self.safety = tk.Scale(self.master, from_=200, to=-200, orient=tk.VERTICAL, length=300, tickinterval=50,
                               bg='#ff474c',
                               font='bold', fg='white')
        self.safety.set(100)
        self.safety.grid(row=3, column=7, sticky=tk.W + tk.E)
        self.over = False

    def get_values(self):
        self.coeff = {'sbio': self.structural_info.get() / 100, 'sdrug': self.structural_drug.get() / 100,
                      'chem': self.chemistry.get() / 100, 'bio': self.biology.get() / 100,
                      'dise': self.disease.get() / 100, 'gen': self.genetic.get() / 100,
                      'info': self.information.get() / 100, 'safe': self.safety.get() / 100}
        self.over = True


class config_gui:
    def __init__(self, get_only=False):
        self.new = False
        self.active = True
        self.config = configparser.ConfigParser()
        self.config_file_path = Path('~/.targetdb/config.ini').expanduser()
        self.config_file_path.parent.mkdir(exist_ok=True, parents=True)

        if self.config_file_path.is_file():
            self.config.read(str(self.config_file_path))
            if not get_only:
                self.config_panel(new=False)
        else:
            self.config['database_path'] = {}
            self.config['output_path'] = {}
            self.config['pubmed_email'] = {}
            self.config['executables'] = {}
            self.config['database_path']['targetdb'] = ''
            self.config['database_path']['chembl'] = ''
            self.config['pubmed_email']['email'] = ''
            self.config['output_path']['db_files'] = ''
            self.config['output_path']['single'] = ''
            self.config['output_path']['list'] = ''
            self.config['executables']['blast'] = ''
            self.config['executables']['blastdb_path'] = ''
            self.config['executables']['fpocket'] = ''
            self.config_panel(new=True)

    def config_panel(self, new=False):
        self.new = new
        self.config_window = tk.Toplevel()
        self.w_header = tk.Label(self.config_window, text='CONFIGURATION', bg="#bbf7b4", font='Helvetica 16 bold')
        self.w_section_database = tk.Label(self.config_window, text='Databases', relief="raised", bg='#e6e6e6')
        self.w_section_outputs = tk.Label(self.config_window, text='Outputs', relief="raised", bg='#e6e6e6')
        self.w_section_pubmed = tk.Label(self.config_window, text='Pubmed', relief="raised", bg='#e6e6e6')
        self.message_pubmed = tk.StringVar()
        self.message_pubmed.set(
            """Enter your email address (used for pubmed searches - pubmed api requires an email for batch requests)""")
        self.w_pubmed_msg = tk.Message(self.config_window, textvariable=self.message_pubmed, justify='center',
                                       width=500)
        self.targetdb_path = tk.StringVar()
        self.w_targetdb = tk.Entry(self.config_window, textvariable=self.targetdb_path, width=100)
        self.list_path = tk.StringVar()
        self.w_list = tk.Entry(self.config_window, textvariable=self.list_path, width=100)
        self.single_path = tk.StringVar()
        self.w_single = tk.Entry(self.config_window, textvariable=self.single_path, width=100)
        self.pubmed_email = tk.StringVar()
        self.w_pubmed_email = tk.Entry(self.config_window, textvariable=self.pubmed_email, width=100)

        self.w_header.grid(row=0, column=0, columnspan=5, sticky="ew")
        self.w_section_database.grid(row=1, column=0, columnspan=2, sticky="ew")
        tk.Message(self.config_window, text="""Please browse to a valid targetDB database file:""", width=500).grid(
            row=2, column=0)
        self.w_targetdb.grid(row=2, column=1, columnspan=2, sticky="ew")
        tk.Button(self.config_window, text="Browse...", bg="#a0a0a0", command=self.select_targetdb).grid(row=2,
                                                                                                         column=4)

        self.w_section_outputs.grid(row=3, column=0, columnspan=2, sticky="ew")

        tk.Message(self.config_window, text="""Please select a folder in which to save LISTS outputs:""",
                   width=500).grid(row=5, column=0)
        self.w_list.grid(row=5, column=1, columnspan=2, sticky="ew")
        tk.Button(self.config_window, text="Browse...", bg="#a0a0a0", command=self.select_folder_list).grid(row=5,
                                                                                                            column=4)

        tk.Message(self.config_window, text="""Please select a folder in which to save SINGLE outputs:""",
                   width=500).grid(row=7, column=0)
        self.w_single.grid(row=7, column=1, columnspan=2, sticky="ew")
        tk.Button(self.config_window, text="Browse...", bg="#a0a0a0", command=self.select_folder_single).grid(row=7,
                                                                                                              column=4)

        self.w_section_pubmed.grid(row=8, column=0, columnspan=2, sticky="ew")
        self.w_pubmed_msg.grid(row=9, column=0)
        self.w_pubmed_email.grid(row=9, column=1, columnspan=3, sticky="ew")

        tk.Button(self.config_window, text="Save & Close", bg="#00cc99", command=self.save,
                  font='Helvetica 12 bold').grid(row=10, column=0, columnspan=4, sticky="ew")

        if not new:
            self.targetdb_path.set(self.config['database_path']['targetdb'])
            self.list_path.set(self.config['output_path']['list'])
            self.single_path.set(self.config['output_path']['single'])
            self.pubmed_email.set(self.config['pubmed_email']['email'])

        if self.new:
            self.config_window.protocol('WM_DELETE_WINDOW', self.confirm)
        else:
            self.config_window.protocol('WM_DELETE_WINDOW', self.close_window)
        while self.active:
            self.config_window.update()

    def confirm(self):
        if askokcancel('Quit',
                       'No configuration file is present.\nIf you close the window without setting the parameters and save them targetDB will close'):
            self.active = False
            self.config_window.destroy()

    def close_window(self):
        self.config_window.destroy()
        self.active = False

    def save(self):
        self.config['database_path']['targetdb'] = self.w_targetdb.get()
        self.config['output_path']['list'] = self.w_list.get()
        self.config['output_path']['single'] = self.w_single.get()
        self.config['pubmed_email']['email'] = self.w_pubmed_email.get()
        if self.check_config():
            with self.config_file_path.open(mode='w') as cfile:
                self.config.write(cfile)
            self.config_window.destroy()
            self.active = False
        else:
            return False

    def select_targetdb(self):
        targetDB_path = Path(askopenfilename(title='Select TargetDB sqlite database',
                                             filetypes=[("sqliteDB", "*.db")]))
        if targetDB_path.is_file():
            self.targetdb_path.set(targetDB_path)
        else:
            showerror('Error', 'Please select a valid file')

    def select_folder_list(self):
        folder = Path(askdirectory(title='Select directory to save files'))
        if folder.is_dir():
            self.list_path.set(folder)
        else:
            showerror('Error', 'Please select a valid folder')

    def select_folder_single(self):
        folder = Path(askdirectory(title='Select directory to save files'))
        if folder.is_dir():
            self.single_path.set(folder)
        else:
            showerror('Error', 'Please select a valid folder')

    def check_config(self):
        error = False
        msg = ''
        if not Path(self.config['database_path']['targetdb']).is_file():
            error = True
            msg += 'TargetDB database file not found, please select a valid file\n'
        if not Path(self.config['output_path']['list']).is_dir():
            error = True
            msg += 'Folder for LIST output not valid, please select a valid folder\n'
        if not Path(self.config['output_path']['single']).is_dir():
            error = True
            msg += 'Folder for SINGLE output not valid, please select a valid folder\n'
        if not is_email(self.config['pubmed_email']['email']):
            error = True
            msg += 'Email Address not valid, please enter a valid email address\n'
        if error:
            showerror('Error - Invalid parameters', msg)
            return False
        else:
            return True


def is_email(email):
    match = re.match('^[_a-z0-9-]+(\.[_a-z0-9-]+)*@[a-z0-9-]+(\.[a-z0-9-]+)*(\.[a-z]{2,4})$', email)
    if match is None:
        return False
    return True


class targetDB_gui:
    def __init__(self):
        self.window = tk.Tk()
        dr.get_global_param()
        self.targetDB_path = dr.targetDB

        # ================================================================================ #
        # ============================== WIDGET BUILDING ================================= #
        # ================================================================================ #

        self.title = tk.Label(self.window, height=2, text='TargetDB', bg='#00cc99', font='Helvetica 16 bold')

        self.mode = tk.StringVar()

        self.mode.set('list')
        self.mode_label = tk.Label(self.window, text='Mode', relief="raised", bg='#e6e6e6')
        self.mode_single = tk.Radiobutton(self.window, text="Single mode", variable=self.mode, value='single')
        self.mode_list = tk.Radiobutton(self.window, text="List mode", variable=self.mode, value='list')
        self.mode_spider = tk.Radiobutton(self.window, text="Spider plot only", variable=self.mode, value='spider')

        self.target_input_label = tk.Label(self.window, text='Targets list', relief="raised", bg='#e6e6e6')
        self.message_target = tk.StringVar()
        self.message_target.set(
            """Please enter the list of targets. You can either enter a list of gene names directly from excel or a list of comma separated values (eg. "DYRK1A,CYP3D6")""")
        self.target_message = tk.Message(self.window, textvariable=self.message_target, justify="center")
        self.input_text = ScrolledText(self.window)

        self.settings_button = tk.Button(self.window, text='Settings', bg='#b7b7b7', command=self.get_settings,
                                         font='Helvetica 10 bold')
        self.go_button = tk.Button(self.window, text='Start', bg='#00cc99', command=self.launch
                                   , font='Helvetica 12 bold')

        self.close = tk.Button(self.window, text='Exit', fg='red', command=self.window.destroy
                               , font='Helvetica 12 bold')

        # ================================================================================ #
        # =========================== GRID CONSTRUCTION ================================== #
        # ================================================================================ #

        self.title.grid(row=0, column=0, columnspan=4, sticky="ew")
        self.mode_label.grid(row=1, column=0, ipadx=10, pady=10, sticky="w")
        self.mode_single.grid(row=1, column=1, pady=10)
        self.mode_list.grid(row=1, column=2, pady=10)
        self.mode_spider.grid(row=1, column=3, pady=10)

        self.target_input_label.grid(row=3, column=0, ipadx=10, sticky="wne")
        self.target_message.grid(row=4, column=0, rowspan=2, sticky='nw')
        self.settings_button.grid(row=6, column=0, sticky='nwe', ipadx=10)
        self.input_text.grid(row=3, column=1, columnspan=3, rowspan=4)

        self.go_button.grid(row=8, column=1, columnspan=2, sticky="ew", pady=10, padx=10)
        self.close.grid(row=8, column=3, sticky='ew', ipadx=20, pady=10, padx=10)

        # ================================================================================ #
        # =========================== LAUNCHING THE APP ================================== #
        # ================================================================================ #
        self.window.mainloop()

    def get_settings(self):
        self.settings = config_gui(get_only=False)

    def launch(self):
        target_list = self.input_text.get("1.0", tk.END).splitlines()
        mode = self.mode.get()
        if len(target_list) == 1:
            target_list = target_list[0].split(',')
        if len(target_list) == 1:
            showerror('Error'
                      ,
                      'Your list of targets is empty or does not seems to be copy/paste from excel or a comma separated string (eg. DYRK1A,PREP,BACE1,...) ')
            return None
        self.gene_df = g2id.gene_to_id(target_list, targetDB_path=self.targetDB_path)
        message = ''
        if mode == 'list':
            dr.get_list_excel(self.gene_df)
        if mode == 'single':
            for gene_name in self.gene_df.index:
                message += dr.get_single_excel(self.gene_df.loc[gene_name])
            info_message(message)
        if mode == 'spider':
            self.top_spider = tk.Toplevel(self.window)
            self.top_spider.protocol("WM_DELETE_WINDOW", self.closing_top)
            self.top_spider.geometry("600x800")
            self.top_spider.resizable(False, True)

            spider_window = ScrolledFrame(self.top_spider)
            spider_window.pack(expand=True, fill='both')
            # grid_width = 3
            # col = 0
            # row = 0
            for gene_name in self.gene_df.index:
                target_desc = td.get_descriptors_list(self.gene_df.loc[gene_name]['uniprot_ids'][0],
                                                      targetdb=self.targetDB_path)
                tscore = td.target_scores(target_desc, mode='single')
                target_desc = target_desc.merge(tscore.scores, on='Target_id', how='left')
                score_col = ['structure_info_score', 'structural_drug_score', 'chemistry_score', 'biology_score',
                             'disease_score', 'genetic_score', 'information_score', 'safety_score']
                target_score = target_desc[score_col] * 10
                target_score.index = target_desc.Target_id
                target_score.fillna(0, inplace=True)
                target_score = target_score.rename(columns={'structure_info_score': 'Structural information',
                                                            'structural_drug_score': 'Structural Druggability',
                                                            'chemistry_score': 'Chemistry', 'biology_score': 'Biology',
                                                            'disease_score': 'Diseases',
                                                            'genetic_score': 'Genetic Association',
                                                            'information_score': 'Literature',
                                                            'safety_score': 'Safety'})
                self.spider_plot = self.make_spider_plot \
                    (target_score.loc[self.gene_df.loc[gene_name]['uniprot_ids'][0]].values, target_score.columns,
                     target_name=self.gene_df.loc[gene_name]['symbol'])
                canvas = FigureCanvasTkAgg(self.spider_plot, master=spider_window.inner)
                canvas.get_tk_widget().pack(expand=True, fill='both')
                # canvas.get_tk_widget().grid(column=col,row=row,sticky='ew')
                # col+=1
                # if col == grid_width:
                #     col = 0
                #     row+=1
                plt.close(self.spider_plot)
            self.window.wait_window(self.top_spider)

    def closing_top(self):
        self.top_spider.destroy()

    def make_spider_plot(self, data, labels, target_name=''):
        fig = plt.figure(figsize=(4, 3))
        ax = fig.add_axes([0, 0, 0.6, 0.8], polar=True)
        ax.spines['polar'].set_visible(False)
        N = len(data)
        bar_width = 2 * np.pi / N
        theta = np.arange(0, 2 * np.pi, 2 * np.pi / N)
        colors = ['#95d0fc', '#0485d1', '#b790d4', '#87ae73', '#fec615', '#fb7d07', '#95a3a6', '#ccff33']
        ax.bar(0, 1, bottom=9, width=2 * np.pi, color='r', linewidth=0, alpha=0.3)
        ax.bar(0, 5, bottom=4, width=2 * np.pi, color='lime', linewidth=0, alpha=0.2)
        ax.bar(0, 3, bottom=1, width=2 * np.pi, color='gold', linewidth=0, alpha=0.2)
        ax.bar(0, 1, width=2 * np.pi, color='r', linewidth=0)
        for i in range(len(data)):
            ax.bar(theta[i], data[i], width=bar_width, align='center', label=list(labels)[i], color=colors[i],
                   edgecolor='black', linewidth=1.5)
        plt.title(target_name, weight='bold', fontsize=14)
        ax.set_xticks(theta)
        x_labels = [''] * len(theta)
        ax.set_xticklabels(x_labels)
        ax.yaxis.grid(True)
        ax.xaxis.grid(False)
        fig.legend(loc=7, fontsize=10, fancybox=True, markerscale=1.2)
        ax.set_yticks([])
        return fig


def error(message):
    showerror("Error", message)


def info_message(message):
    showinfo("Message", message)
