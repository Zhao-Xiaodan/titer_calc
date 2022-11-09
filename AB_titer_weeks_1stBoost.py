import os
import pathlib
from operator import itemgetter, attrgetter
import datetime

# from cv2 import cv2

import numpy as np
from math import sqrt
import matplotlib.pyplot as plt
import seaborn as sns
import pandas as pd
from scipy import interpolate

import operator

def get_truth(week_get, relate, week_target):
    return relate(week_get, week_target)

def SelectPpPairedWeek(pairs_ls, dataFrame_1stBoost):
    for pairs in pairs_ls:
        pp, relate, week_target = pairs
        condition_0 = get_truth(dataFrame_1stBoost['Days from 1st Booster'],
                                relate, week_target)
        condition = condition_0 & (dataFrame_1stBoost['Participant'].isin([pp]))

        #  print(dataFrame_1stBoost['Participant'].isin([pp]))
        #  print(sum(dataFrame_1stBoost['Participant'].isin([pp])))
        #  breakpoint()
        dataFrame_1stBoost =  dataFrame_1stBoost.loc[~condition]
    return dataFrame_1stBoost

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
    plt.ylim(-1500, 26500)
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
    plt.ylim(-1500, 26500)

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
    sns.set_theme(style='darkgrid')
    #  dataFitting = dataFrame[[\
           #  'Participant', byColumn, 'Titer']].sort_values(
            #  by=['Participant', byColumn])
    #  dataFitting.to_csv("fvaccinated.csv", index=False)

    #  #  customPalette = sns.color_palette("Set1", as_cmap = True)
    dataFrame = days2weeks(dataFrame, byColumn, 7)
    dataFrame = selectParticipantByCount(dataFrame, 1)
    if False:
        sns.lineplot(
            data=dataFrame,
            #  x="Days from 2nd Jab", y="Titer", hue="vaccineCurrent", style="dataInfection",
            x=byColumn, y="Titer", hue="vaccineCurrent", marker='o',
            markers=True, dashes=False,units="Participant", estimator=None, lw=1,
            markersize=8
        )

    textLabel(dataFrame, byColumn)
    textLabel_begin(dataFrame, byColumn)
    BSpline(dataFrame, byColumn)
    plt.ylim(-3000, 34000)
    plt.show()
    #


def BSpline(dataFrame, byColumn):
    dataFrame = dataFrame.loc[dataFrame[byColumn] != -np.inf]
    dataFrame = dataFrame.sort_values(by=['Participant', byColumn])
    dataFrame = dataFrame[ ['Participant', byColumn, 'Titer', 'vaccineCurrent'] ]


    gb = dataFrame[['Participant', byColumn, 'Titer', 'vaccineCurrent']].groupby('Participant')

    fig, ax = plt.subplots(figsize=(5,4))
    #  plt.figure()
    #  plt.clf()
    #  ax = plt.gca()
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


# ============ main ===============

df = pd.read_csv('titerMaster.csv')
#  df = pd.read_excel('readTiterMaster_06092022.xlsx', header=0, engine='openpyxl')

# cap the max titer to be 25000
df.loc[df['Titer']>25000, 'Titer'] = 25000

# only analyse 'WT RBD'
df= df[df['Varian'] == 'WT RBD']
df = df.loc[df["vaccineCategory"] == "mRNA"]

ppSelection_Pfizer = [32, 121, 29]
#  ppSelection_Pfizer = [10, 18, 29, 32, 95, 121]
ppSelection_Moderna = [9, 14, 23, 70, 89, 138,  143, 144, 174] + [68, 120, 103]
#  ppSelection_Moderna = [9, 14, 23, 70, 89, 138, 143, 144, 174] + [7, 68, 120, 169, 103, 39, 65, 113]
#  ppSelection_Moderna = [9, 14, 23, 70, 89, 138, 143, 144, 174] + [7, 68, 120, 169, 103]
ppSelection = ppSelection_Pfizer + ppSelection_Moderna
df = df.loc[~df['Participant'].isin(ppSelection)]


# select pp whose first Test Date is within certain weeks only select fVaccinated data
ind_min = df[['Participant', 'Days from 2nd Jab']].groupby('Participant').idxmin()
ind_min = ind_min.values.flatten()
df_temp = df[['Test Date', 'Participant', 'Days from 2nd Jab']]
df_temp.loc[ind_min, 'idxmin_mask'] = 1
#  df_temp.to_csv("temp_freq_count.csv")
df_temp.loc[(df_temp['idxmin_mask'] == 1) & (df_temp['Days from 2nd Jab'] < 5*7), 'idxmin_mask2'] = 1
ind_pp = df_temp.loc[df_temp['idxmin_mask2'] ==1]['Participant'].values
#  df_temp.to_csv("temp_freq_count.csv")
#  df = df.loc[df['Participant'].isin(ind_pp)] # comment off for 1st Booster data
#  df.to_csv("df_pp.csv")

#  breakpoint()
#  vaccineList = ['Moderna']
#  vaccineList = ['Pfizer']
vaccineList = ['Moderna', 'Pfizer']
#  multiLinePlot(vaccineList[0])

# mask infected
df['dataInfection'] = np.where(df['daysInfection']>0, 1, 0)
# mask fully vaccinated
df.loc[df["Days from 1st Booster"]<0, "fVaccinated"] = 1
# mask 1st Booster
df.loc[df['Days from 1st Booster'] != -np.inf, 'data1stBooster'] = 1
df.loc[df['Days from 2nd Booster'] >= 0, 'data1stBooster'] = 0

# ==== select 1st Booster only include one negative value if it has ====

#  ind1stBooster_max = df.loc[(df["Days from 1st Booster"] != -np.inf) & \
#                         #  (df["daysInfection"] == -np.inf) &
#                         (df["Days from 1st Booster"]<0)]\
#      [["Participant", "Days from 1st Booster"]].groupby("Participant").idxmax()
#  ind1stBooster_max = ind1stBooster_max.values.flatten()
#  df.loc[ind1stBooster_max, "data1stBooster"] = 1
#  #  df.loc[((df["Days from 1st Booster"]>0) & (df["Days from 2nd Booster"]<0) & \
#          #  (df["daysInfection"]<0)), "data1stBooster"] = 1
#  df.loc[((df["Days from 1st Booster"]>0) & (df["Days from 2nd Booster"]<0)),\
#           "data1stBooster"] = 1

# mask 2nd Booster
ind2ndBooster = df.loc[(df["Days from 2nd Booster"] != -np.inf) & \
                       (df["daysInfection"] == -np.inf) &
                       (df["Days from 2nd Booster"]<0)]\
    [["Participant", "Days from 1st Booster"]].groupby("Participant").idxmax()
ind2ndBooster = ind2ndBooster.values.flatten()
df.loc[ind2ndBooster, "data2ndBooster"] = 1
#  df.loc[(df["Days from 2nd Booster"]>0) & (df["daysInfection"]<0), "data2ndBooster"] = 1
df.loc[df["Days from 2nd Booster"]>0 , "data2ndBooster"] = 1

# select the latest data from the duplicated ones
df = df.sort_values(by=['Participant', 'Test Date', 'Varian', 'File Number'])
df = df.groupby(['Participant', 'Test Date', 'Varian']).tail(1).reset_index(drop=True)

# export raw data in wide format
df_wide = df.pivot(columns='Participant',
              values=['Days from 2nd Jab', 'Titer']
              ).swaplevel(0, 1, axis=1).sort_index(1)
df_wide = df_wide.apply(lambda x: pd.Series(x.dropna().values))
df_wide.to_csv("export_1stBoost.csv")
# =======================================================================

# partition data based on above masking
dataFrame_fVaccinated, dataFrame_1stBoost, \
    dataFrame_2ndBoost = \
    mRNAvaccine(df)
# define range of data to be viewed
dataFrame_fVaccinated = dataFrame_fVaccinated.loc[dataFrame_fVaccinated[
'Days from 2nd Jab']<200]
dataFrame_1stBoost = dataFrame_1stBoost.loc[
    (dataFrame_1stBoost[
    'Days from 1st Booster']<200) &
    (dataFrame_1stBoost['Days from 1st Booster'] > -100)]


pairs_ls = [
    (15, operator.lt, -9),
    (112, operator.lt, -9),
    (49, operator.lt, -9),
    (18, operator.lt, -7.5),
    (119, operator.lt, -7.5),
    (7, operator.eq, -8),
    (59, operator.gt, 35),
]

dataFrame_1stBoost = SelectPpPairedWeek(pairs_ls, dataFrame_1stBoost)

allPP = dataFrame_1stBoost['Participant'].unique()
dataFrame_1stBoost.loc[~dataFrame_1stBoost['Participant'].isin(allPP)]
#  breakpoint()
sns.set_theme(style='darkgrid')
#  sns.set_theme(style='whitegrid')
#  oneClickPlot(dataFrame_fVaccinated, 'Days from 2nd Jab')
oneClickPlot(dataFrame_1stBoost, 'Days from 1st Booster')

breakpoint()
list(map(onePlotAll, vaccineList))
#  list(map(onePlot, vaccineList))
breakpoint()

