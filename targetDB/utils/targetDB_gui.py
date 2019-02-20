#!/usr/bin/env python
import tkinter as tk
from tkinter.scrolledtext import ScrolledText
from tkinter.messagebox import showerror
import matplotlib.pyplot as plt
from matplotlib.backends.backend_tkagg import FigureCanvasTkAgg
import numpy as np

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
        self.message_target.set("""Please enter the list of targets. You can either enter a list of gene names directly from excel or a list of comma separated values (eg. "DYRK1A,CYP3D6")""")
        self.target_message = tk.Message(self.window, textvariable=self.message_target, justify="center")
        self.input_text = ScrolledText(self.window)

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
        self.input_text.grid(row=3, column=1, columnspan=3, rowspan=4)

        self.go_button.grid(row=8, column=1, columnspan=2, sticky="ew", pady=10, padx=10)
        self.close.grid(row=8, column=3, sticky='ew', ipadx=20, pady=10, padx=10)

        # ================================================================================ #
        # =========================== LAUNCHING THE APP ================================== #
        # ================================================================================ #
        self.window.mainloop()

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

        if mode == 'list':
            dr.get_list_excel(self.gene_df)
        if mode == 'single':
            for gene_name in self.gene_df.index:
                dr.get_single_excel(self.gene_df.loc[gene_name])
        if mode == 'spider':
            self.top_spider = tk.Toplevel(self.window)
            self.top_spider.protocol("WM_DELETE_WINDOW", self.closing_top)
            self.top_spider.geometry("600x800")
            self.top_spider.resizable(False, True)

            spider_window = ScrolledFrame(self.top_spider)
            spider_window.pack(expand=True, fill='both')

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
