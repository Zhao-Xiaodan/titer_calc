#  import os
#  import pathlib
#  from operator import itemgetter, attrgetter
#  import datetime
#
#  from cv2 import cv2
#
import statannot
import scipy
import numpy as np
#  from math import sqrt
import matplotlib.pyplot as plt
import matplotlib as mpl
import seaborn as sns
import pandas as pd
from scipy import interpolate

def convert_pvalue_to_asterisks(pvalue):
    if pvalue <= 0.0001:
        return "****"
    elif pvalue <= 0.001:
        return "***"
    elif pvalue <= 0.01:
        return "**"
    elif pvalue <= 0.05:
        return "*"
    return "ns"

def selectParticipantByCount(df, frequence):
    counts = df['Participant'].value_counts()
    df = df[df['Participant'].isin(counts[counts > frequence].index)]
    # an alternative way below
    #  df['freq_count'] = df.groupby('Participant')['Participant'].transform('count')
    #  df = df.loc[df['freq_count'] > 1]
    return df

def skipMonthWindow(df, monthWindow=[-2, 1, 4]):
    df.loc[:, 'month after 2nd Jab'] = df['Days from 2nd Jab'].floordiv(30)

    #  df['month after 2nd Jab'] = df['Days from 2nd Jab'].floordiv(30)
    df = df.loc[~df["month after 2nd Jab"].isin(monthWindow)]

def groupPlot(dataFrame, dataFrame_1stBoost, dataFrame_2ndBoost, \
              dataFrame_infection, xData):
    fig, ax = plt.subplots()
    dataFrame_infection.groupby("Participant").plot(x=xData, y="Titer", \
                                                    marker="o", color='k', markersize=15, \
                                                    linewidth=0, ax=ax)
    dataFrame.groupby("Participant").plot(x=xData, y="Titer", \
                                          marker="o", markersize=10, \
                                          linewidth=0.3, ax=ax)
    dataFrame_1stBoost.groupby("Participant").plot(x=xData, y="Titer", \
                                                   marker="*", color='w', linewidth=0,\
                                                   markersize=6, ax=ax)

    dataFrame_2ndBoost.groupby("Participant").plot(x=xData, y="Titer", \
                                                   marker="^", color='w', linewidth=0,\
                                                   markersize=15, ax=ax)
    #  breakpoint()
    #  plt.ylim(-1500, 26500)
    ax.get_legend().remove()
    plt.ylabel('Titer - times of dilution')


def mRNAvaccine(dataFrame):
    #  vaccineCategory
    dataFrame_fVaccinated = dataFrame.loc[dataFrame['fVaccinated'] == 1]
    dataFrame_1stBoost = dataFrame.loc[dataFrame['data1stBooster'] == 1]
    dataFrame_2ndBoost = dataFrame.loc[dataFrame['data2ndBooster'] == 1]
    return dataFrame_fVaccinated, dataFrame_1stBoost, dataFrame_2ndBoost

def vlPlot(dataFrame):
    data_roi = dataFrame[["Days from 2nd Jab", "Titer"]]
    pivoted = data_roi.pivot(columns="Days from 2nd Jab", values="Titer")
    sns.violinplot(data=pivoted, palette="Set3", bw=.2, cut=1, linewidth=1)

def linePlot(dataFrame):
    sns.lineplot(x="Days from 2nd Jab", y="Titer",
                 data=dataFrame)
    #  plt.ylim(-1500, 26500)

def days2weeks(dataFrame, byColumn, period):
    dataFrame.loc[:,byColumn] = dataFrame[byColumn].floordiv(period)
    return dataFrame

def textLabel_begin(dataFrame, byColumn):
    Txt=dataFrame[ ['Participant', \
                    byColumn, \
                    'Titer'] ]

    indTxt=Txt.groupby('Participant')[byColumn].idxmin()
    indTxt=indTxt.values.flatten()
    Txt = Txt.loc[indTxt].values
    for i in range(len(Txt)):
        plt.text(Txt[i][1], Txt[i][2], str(int(Txt[i][0])))

def textLabel(dataFrame, byColumn):
    Txt=dataFrame[ ['Participant', \
                    byColumn, \
                    'Titer'] ]

    indTxt=Txt.groupby('Participant')[byColumn].idxmax()
    indTxt=indTxt.values.flatten()
    Txt = Txt.loc[indTxt].values
    for i in range(len(Txt)):
        plt.text(Txt[i][1], Txt[i][2], str(int(Txt[i][0])))

def oneClickPlot(dataFrame, byColumn):

    dataFitting = dataFrame[[\
                             'Participant', byColumn, 'Titer']].sort_values(
                                 by=['Participant', byColumn])
    dataFitting.to_csv("fVaccinated.csv", index=False)

    #  #  customPalette = sns.color_palette("Set1", as_cmap = True)
    dataFrame = days2weeks(dataFrame, byColumn, 7)
    dataFrame = selectParticipantByCount(dataFrame, 1)
    sns.lineplot(
        data=dataFrame,
        #  x="Days from 2nd Jab", y="Titer", hue="vaccineCurrent", style="dataInfection",
        x=byColumn, y="Titer", hue="vaccineCurrent", marker='o',
        markers=True, dashes=False,units="Participant", estimator=None, lw=0,
        markersize=8
    )

    #  textLabel(dataFrame, byColumn)
    #  textLabel_begin(dataFrame, byColumn)
    BSpline(dataFrame, byColumn)
    plt.show()
    #


def BSpline(dataFrame, byColumn):
    dataFrame = dataFrame.loc[dataFrame[byColumn] != -np.inf]
    dataFrame = dataFrame.sort_values(by=['Participant', byColumn])
    dataFrame = dataFrame[ ['Participant', byColumn, 'Titer', 'vaccineCurrent'] ]


    gb = dataFrame[['Participant', byColumn, 'Titer', 'vaccineCurrent']].groupby('Participant')

    plt.figure()
    plt.clf()
    ax = plt.gca()
    for Participant, frame in gb:
        numData = frame['Participant'].count()
        if numData == 2: kv=1
        elif numData > 2: kv=2
        #  elif numData > 3: kv=3

        if numData > 1: # not apply for single data point
            color = next(ax._get_lines.prop_cycler)['color']
            #  print(numData, kv)
            x = frame[byColumn].values.tolist()
            y = frame['Titer'].values.tolist()

            tck = interpolate.splrep(x, y, s=0, k=kv)
            x_new = np.linspace(min(x), max(x), 100)
            y_fit = interpolate.BSpline(*tck)(x_new)
            if pd.isnull(y_fit).any():
                print(frame['Participant'])
            #  plt.plot(x, y, 'o', color=color)
            #  plt.plot(x_new, y_fit, '-', linewidth=1, color=color)
            if frame['vaccineCurrent'].str.contains('Pfizer').any():
                plt.plot(x, y, 'o', color=color)
                plt.plot(x_new, y_fit, '-', linewidth=1, color=color)
                #  plt.plot(x_new, y_fit, '-', linewidth=1, color='C0')
            elif frame['vaccineCurrent'].str.contains('Moderna').any():
                plt.plot(x, y, '*', markersize=10, color=color)
                plt.plot(x_new, y_fit, '-', linewidth=1, color=color)
                #  plt.plot(x_new, y_fit, '-', linewidth=1, color='C1')

def removeDataPointbyList_Xrange(pp_list, dataFrame, byColumn, Xrange):
    condition = (dataFrame[byColumn] < Xrange) & \
        (dataFrame_1stBoost['Participant'].isin(pp_list))
    dataFrame =  dataFrame.loc[~condition]
    return dataFrame

# ============ main ===============

#  df = pd.read_excel('readTiterMaster_20092022.xlsx', header=0, engine='openpyxl')

df = pd.read_csv('titerMaster.csv')
# cap the max titer to be 25000
df.loc[df['Titer']>25000, 'Titer'] = 25000

# select by same date and same pp
# the selection is not perfect, need t consider later
fileList = df.loc[(df['Varian'] == 'BA4/5 RBD') |
                  (df['Varian'] == 'BA1/2 RBD')]['File Number'].values
#  df = df.loc[df['File Number'].isin(fileList)]
# instead of trim df, just add a mask column then can retrieve back other data
df.loc[df['File Number'].isin(fileList), 'varianTest'] = 1
#  pp = df.loc[df['Varian'] == 'BA4/5 RBD']['Participant'].values
#  df = df.loc[(df['Test Date'].isin(testDate)) & (df['Participant'].isin(pp))]
df = df.sort_values(by=['Test Date', 'Participant'])

# mask infected
df['dataInfection'] = np.where(df['daysInfection']>0, 1, 0)
# mask fully vaccinated
df.loc[df["Days from 1st Booster"]<0, "fVaccinated"] = 1
# mask 1st Booster
df.loc[(df['Days from 1st Booster'] >= 0) & \
       (df['Days from 2nd Booster'] < 0), 'data1stBooster'] = 1
# mask 2nd Booster
df.loc[df['Days from 2nd Booster'] >= 0, 'data2ndBooster'] = 1

# Zhou Yu select participants
selection2dose = [18, 19, 32, 39, 70, 15, 117, 35, 73, 4]
#  selection2dose = sorted(selection2dose)
selection2doseInfection = [262, 263, 95]
#  selection2doseInfection= sorted(selection2doseInfection)
selection3dose = [120, 119, 32, 70, 15, 161, 127, 103, 23, 162]
#  selection3dose = sorted(selection3dose)
selection3doseInfection = [9, 14, 70, 255, 19, 174, 266, 261, 32, 118, 267]
#  selection3doseInfection = sorted(selection3doseInfection)
selection4dose = [121, 138, 259, 23, 143, 144, 271, 62]
#  selection4dose = sorted(selection4dose)
selection3sinoInfection = [177, 256, 94, 265]
#  selection3sinoInfection = sorted(selection3sinoInfection)

df_2dose = pd.DataFrame({'2dose':selection2dose})
df_2doseInfe = pd.DataFrame({'2doseInfe':selection2doseInfection})
df_3dose = pd.DataFrame({'3dose':selection3dose})
df_3doseInfe = pd.DataFrame({'3doseInfe':selection3doseInfection})
df_4dose = pd.DataFrame({'4dose':selection4dose})
df_3sinoInfe = pd.DataFrame({'3sinoInfe':selection3sinoInfection})

df_selection = pd.concat(
    [df_2dose, df_2doseInfe, df_3dose, df_3doseInfe,
           df_3sinoInfe, df_4dose], axis=1)
for col in df_selection:
    df_selection[col] = df_selection[col].sort_values(ignore_index=True)

# mask the Vaccine Infection History
df.loc[(df['dataInfection'] == 0) &
       (df['Participant'].isin(selection2dose)) &
       (df['Varian'] == 'WT RBD') &
       #  (df['Days from 2nd Jab'] >= 14) &
       #  (df['Days from 2nd Jab'] < 42) &
       (df['fVaccinated'] == 1),
       'Vaccine Infection History'] = '2 dose'

df.loc[(df['dataInfection'] == 0) &
       #  (df['varianTest'] == 1) &
       (df['Participant'].isin(selection2dose)) &
       #  (df['Days from 2nd Jab'] >= 14) &
       #  (df['Days from 2nd Jab'] < 42) &
       (df['fVaccinated'] == 1),
       'Vaccine Infection History'] = '2 dose'

df.loc[(df['dataInfection'] == 1) &
       #  (df['varianTest'] == 1) &
       (df['Participant'].isin(selection2doseInfection)) &
       #  (df['daysInfection'] >= 14) &
       #  (df['daysInfection'] < 42) &
       (df['fVaccinated'] == 1),
       'Vaccine Infection History'] = '2 dose + Infection'

df.loc[(df['dataInfection'] == 0) &
       (df['Participant'].isin(selection3dose)) &
       (df['Varian'] == 'WT RBD') &
       #  (df['Days from 1st Booster'] >= 14) &
       #  (df['Days from 1st Booster'] < 42) &
       (df['data1stBooster'] == 1),
       'Vaccine Infection History'] = '3 dose'

df.loc[(df['dataInfection'] == 0) &
       #  (df['varianTest'] == 1) &
       (df['Participant'].isin(selection3dose)) &
       #  (df['Days from 1st Booster'] >= 14) &
       #  (df['Days from 1st Booster'] < 42) &
       (df['data1stBooster'] == 1),
       'Vaccine Infection History'] = '3 dose'

df.loc[(df['dataInfection'] == 1) &
       #  (df['varianTest'] == 1) &
       (df['Participant'].isin(selection3doseInfection)) &
       #  (df['daysInfection'] >= 14) &
       #  (df['daysInfection'] < 42) &
       (df['data1stBooster'] == 1),
       'Vaccine Infection History'] = '3 dose + Infection'

df.loc[(df['dataInfection'] == 0) &
       #  (df['varianTest'] == 1) &
       (df['Participant'].isin(selection4dose)) &
       #  (df['Days from 2nd Booster'] >= 14) &
       #  (df['Days from 2nd Booster'] < 42) &
       (df['data2ndBooster'] == 1),
       'Vaccine Infection History'] = '4 dose'

# Rename the mask when it is mixVaccine
df.loc[(df['vaccineCategory'] == 'mixVaccine') &
       #  (df['varianTest'] == 1) &
       (df['dataInfection'] == 1) &
       #  (df['daysInfection'] >= 14) &
       #  (df['daysInfection'] < 42) &
       (df['Participant'].isin(selection3sinoInfection)),
       'Vaccine Infection History'] = '3 sino + Infection'

col = df.pop('Vaccine Infection History')
df.insert(4, col.name, col)
#  df.sort_values(by=['Vaccine Infection History'])
#  df = df.sort_values(by=['Vaccine Infection History', 'Titer', 'Test Date', 'Participant'])
#  df = df.sort_values(by=['Vaccine Infection History', 'Test Date', 'Participant'])
#  df = df.sort_values(by=['Vaccine Infection History', 'Participant', 'Test Date', 'Titer'])

#
#  breakpoint()

#  print(df.info(verbose=True))
#  df.to_csv('Checking low titer.csv')
#  gb = df.groupby(['Participant', 'Test Date', 'Varian'])['File Number']

df = df.groupby(['Participant', 'Vaccine Infection History', 'Varian']).apply(
    lambda x: x.iloc[-1]) .reset_index(drop=True)

df = df.sort_values(by=['Vaccine Infection History', 'Participant', 'Test Date', 'Varian'])
df['Varian'] = pd.Categorical(df['Varian'], ["WT RBD", "BA1/2 RBD", "BA4/5 RBD"])
#  df = df.sort_values(by=['Vaccine Infection History', 'Titer', 'Test Date', 'Participant'])
#  print(df.info(verbose=True))
df.to_csv('Checking low titer_drop.csv')

# checking plot data points with selections
gb = df.groupby(['Vaccine Infection History', 'Participant'])[[
    'Vaccine Infection History', 'Participant']].apply(
    lambda x: x.iloc[-1]).reset_index(drop=True)

gb_t = gb.pivot(columns='Vaccine Infection History',
                values='Participant')

for col_plot, col_selection in zip(gb_t, df_selection):
    ls_plot = gb_t[col_plot].dropna().values
    ls_selection = df_selection[col_selection].dropna().values
    #  print(ls_plot)
    #  print(ls_selection)
    ls_temp = [x for x in ls_selection if x not in ls_plot]
    if len(ls_temp) > 0:
        print(f"{col_plot}: {ls_temp}")
    #  [print(x in gb_t[col_plot]) for x in df_selection[col_selection]]

# -------- boxplot ----------
if False:
    order = ['2 dose', '2 dose + Infection', '3 dose', '3 dose + Infection', '4 dose', '3 sino + Infection']
    hue_order = ['WT RBD', 'BA1/2 RBD', 'BA4/5 RBD']
    #  sns.set_theme(style='white')
    #  fig, ax = plt.subplots()
    sns.set_theme(style='darkgrid')
    # ------- boxplot ---------

    ax = sns.boxplot(data=df, x="Vaccine Infection History", y="Titer",
                     hue="Varian", order=order, hue_order=hue_order)

    for patch in ax.artists:
        fc = patch.get_facecolor()
        patch.set_facecolor(mpl.colors.to_rgba(fc, 0.4))
    sns.stripplot(data=df, x="Vaccine Infection History", y="Titer",
                  hue="Varian", order=order, hue_order=hue_order, dodge=True, ax=ax)

    # ------------ legend ----------------
    handles, labels = ax.get_legend_handles_labels()
    l = plt.legend(handles[0:4], labels[0:3])
    #  #  l = plt.legend(handles[0:2], labels[0:2], bbox_to_anchor=(1.05, 1), loc=2, borderaxespad=0.)

    plt.yscale('log')
    mng = plt.get_current_fig_manager()
    mng.resize(1680, 1050)
    plt.savefig("fig_varian_comp.png")
    plt.show()

# ---------box plot for fun ---------
if False:
#  if True:
    order = ['2 dose', '2 dose + Infection', '3 dose', '3 dose + Infection', '4 dose', '3 sino + Infection']
    hue_order = ['WT RBD', 'BA1/2 RBD', 'BA4/5 RBD']
    #  sns.set_theme(style='white')
    #  fig, ax = plt.subplots()
    sns.set_theme(style='darkgrid')
    # ------- boxplot ---------
    PROPS = {
    'boxprops':{'facecolor':'none', 'edgecolor':'w', 'lw':2},
    'medianprops':{'color':'w'},
    'whiskerprops':{'color':'w'},
    'capprops':{'color':'w'},
    }

    flierprops = dict(marker='D', markerfacecolor='w', markersize=10,
                  linestyle='none', markeredgecolor='w')

    ax = sns.boxplot(data=df, x="Vaccine Infection History", y="Titer",
                     hue="Varian", order=order, hue_order=hue_order,
                     **PROPS, flierprops=flierprops)

    for patch in ax.artists:
        fc = patch.get_facecolor()
        patch.set_facecolor(mpl.colors.to_rgba(fc, 0.0))
    sns.stripplot(data=df, x="Vaccine Infection History", y="Titer",
                  hue="Varian", order=order, hue_order=hue_order, dodge=True, ax=ax,
                  s=10, marker="D", color='w', linewidth=1, alpha=.7)

    plt.legend([],[], frameon=False)
    #  ax.legend.remove()
    plt.yscale('log')
    mng = plt.get_current_fig_manager()
    mng.resize(1680, 1050)
    plt.savefig("fig_varian_comp.png")
    plt.show()

# -------- ttest and p-value for non-Infection --------
if False:
    dose_Ls = ['2 dose', '3 dose', '4 dose']
    varian_Ls = ['BA1/2 RBD', 'BA4/5 RBD']

    pvalue_dose = []
    for dose_l in dose_Ls:
        stat, pvalue = scipy.stats.ttest_ind(
            (df.loc[(df['Vaccine Infection History'] == dose_l) &
                    (df['Varian'] == 'BA1/2 RBD')]['Titer']
             ),
            (df.loc[(df['Vaccine Infection History'] == dose_l) &
                    (df['Varian'] == 'BA4/5 RBD')]['Titer']
             )
        )
        pvalue_dose.append(pvalue)
        print(pvalue_dose)
    #
    sns.set_theme(style='white')
    fig, ax = plt.subplots()

    #  sns.set_theme(style='darkgrid')
    order = ['2 dose', '3 dose', '4 dose']
    hue_order = ['WT RBD', 'BA1/2 RBD', 'BA4/5 RBD']

    ax = sns.boxplot(data=df, x="Vaccine Infection History", y="Titer",
                     hue="Varian", order=order, hue_order=hue_order)
    for patch in ax.artists:
        fc = patch.get_facecolor()
        patch.set_facecolor(mpl.colors.to_rgba(fc, 0.4))
    sns.stripplot(data=df, x="Vaccine Infection History", y="Titer",
                  hue="Varian", order=order, hue_order=hue_order, dodge=True, ax=ax)

    statannot.add_stat_annotation(
        ax,
        data=df,
        x="Vaccine Infection History",
        y="Titer",
        hue="Varian",
        order=order,
        box_pairs=[
            (("2 dose", "WT RBD"), ("3 dose", "WT RBD")),
            (("2 dose", "WT RBD"), ("4 dose", "WT RBD")),
            (("3 dose", "WT RBD"), ("4 dose", "WT RBD")),
            (("2 dose", "BA1/2 RBD"), ("3 dose", "BA1/2 RBD")),
            (("2 dose", "BA1/2 RBD"), ("4 dose", "BA1/2 RBD")),
            (("3 dose", "BA1/2 RBD"), ("4 dose", "BA1/2 RBD")),
            (("2 dose", "BA4/5 RBD"), ("3 dose", "BA4/5 RBD")),
            (("2 dose", "BA4/5 RBD"), ("4 dose", "BA4/5 RBD")),
            (("3 dose", "BA4/5 RBD"), ("4 dose", "BA4/5 RBD")),
        ],
        test="t-test_ind",
        text_format="star",
        loc="outside",
    )

    handles, labels = ax.get_legend_handles_labels()
    l = plt.legend(handles[0:3], labels[0:3])
    plt.subplots_adjust(left=0.15, right=0.9, top=0.6, bottom=0.1)
    plt.savefig("dose_comapre2.png")
    #  fig.patch.set_visible(False)
    #  ax.axis('off')
    ax.spines['top'].set_visible(False)
    ax.spines['right'].set_visible(False)
    #  ax.spines['bottom'].set_visible(False)
    ax.spines['left'].set_visible(False)
    plt.show()

# -------- ttest and p-value for non-Infection --------
if True:
    dose_Ls = ['2 dose', '3 dose', '4 dose']
    varian_Ls = ['BA1/2 RBD', 'BA4/5 RBD']

    pvalue_dose = []
    for dose_l in dose_Ls:
        stat, pvalue = scipy.stats.ttest_ind(
            (df.loc[(df['Vaccine Infection History'] == dose_l) &
                    (df['Varian'] == 'BA1/2 RBD')]['Titer']
             ),
            (df.loc[(df['Vaccine Infection History'] == dose_l) &
                    (df['Varian'] == 'BA4/5 RBD')]['Titer']
             )
        )
        pvalue_dose.append(pvalue)
        print(pvalue_dose)
    #
    sns.set_theme(style='white')
    fig, ax = plt.subplots()

    #  sns.set_theme(style='darkgrid')
    order = ['2 dose + Infection', '3 dose + Infection', '3 sino + Infection']
    hue_order = ['WT RBD', 'BA1/2 RBD', 'BA4/5 RBD']

    ax = sns.boxplot(data=df, x="Vaccine Infection History", y="Titer",
                     hue="Varian", order=order, hue_order=hue_order)
    for patch in ax.artists:
        fc = patch.get_facecolor()
        patch.set_facecolor(mpl.colors.to_rgba(fc, 0.4))
    sns.stripplot(data=df, x="Vaccine Infection History", y="Titer",
                  hue="Varian", order=order, hue_order=hue_order, dodge=True, ax=ax)

    statannot.add_stat_annotation(
        ax,
        data=df,
        x="Vaccine Infection History",
        y="Titer",
        hue="Varian",
        order=order,
        box_pairs=[
            (("2 dose + Infection", "WT RBD"), ("3 dose + Infection", "WT RBD")),
            (("2 dose + Infection", "WT RBD"), ("3 sino + Infection", "WT RBD")),
            (("3 dose + Infection", "WT RBD"), ("3 sino + Infection", "WT RBD")),
            (("2 dose + Infection", "BA1/2 RBD"), ("3 dose + Infection", "BA1/2 RBD")),
            (("2 dose + Infection", "BA1/2 RBD"), ("3 sino + Infection", "BA1/2 RBD")),
            (("3 dose + Infection", "BA1/2 RBD"), ("3 sino + Infection", "BA1/2 RBD")),
            (("2 dose + Infection", "BA4/5 RBD"), ("3 dose + Infection", "BA4/5 RBD")),
            (("2 dose + Infection", "BA4/5 RBD"), ("3 sino + Infection", "BA4/5 RBD")),
            (("3 dose + Infection", "BA4/5 RBD"), ("3 sino + Infection", "BA4/5 RBD")),
        ],
        test="t-test_ind",
        text_format="star",
        loc="outside",
    )

    handles, labels = ax.get_legend_handles_labels()
    l = plt.legend(handles[0:3], labels[0:3])
    plt.subplots_adjust(left=0.15, right=0.9, top=0.6, bottom=0.1)
    plt.savefig("dose_comapre2.png")
    #  fig.patch.set_visible(False)
    #  ax.axis('off')
    ax.spines['top'].set_visible(False)
    ax.spines['right'].set_visible(False)
    #  ax.spines['bottom'].set_visible(False)
    ax.spines['left'].set_visible(False)
    plt.yscale('log')
    plt.show()

breakpoint()
