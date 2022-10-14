
# Version 2 adds the calculation of CR3022 spike in EC50

import matplotlib.pyplot as plt
import seaborn as sns
import numpy as np
import pandas as pd
from os import listdir, path
from pathlib import Path

# from operator import itemgetter, attrgetter
from datetime import datetime

# load information of participants in excel format
pa = pd.read_excel(
    r"/Users/xiaodanzhao/myProject/AB Data analysis/Participants.xlsx",
    header=1,
    engine="openpyxl",
)

pa = pa.loc[:, pa.notna().any(axis=0)]
# Label rule: (1st and 2nd Jab) - fully vaccinated + 3rd Jab - 1st boost + 4th Jab -2nd boost
jabColNames = ["1st Jab", "2nd Jab", "3rd Jab", "4th Jab"]

# load raw files
rawData_ls = []  # file names of raw data in txt
infoData_ls = []  # file names of info in excel
basepath = Path(__file__).parent
basepath = basepath / "expData"

for i in listdir(basepath):
    if i.endswith(".txt"):
        rawData_ls.append(i)
    elif i.endswith(".xlsx"):
        infoData_ls.append(i)

# seems the sorting is by time, and it is reflected by file number
rawData_ls.sort()
infoData_ls.sort()

data = (
    []
)  # participant no., positive, negative ~ minimum signal of the test, half maximum,  bloodx200, bloodx1000, bloodx5000, bloodx25000, salivax5, salivax20} for each test
data_tag = np.array(["200X", "1000X", "5000X", "25000X", "125000X", "5X", "20X"])
blood_sample_index = [0, 1, 2, 3, 4]
saliva_sample_index = [5, 6]
blood_sample_tag = data_tag[
    blood_sample_index
]  # to extract blood sample dilution factor
# standard_tag = np.array(['400 pM','250 pM','100 pM','50 pM','25 pM','10 pM', '0 pM'])
titerMaster = (
    []
)  # the master file contains date, participant no.,titer no and all relevant test information

noOfRawFiles = len(rawData_ls)


for i in range(noOfRawFiles):
    fileName = rawData_ls[i].strip(".txt")
    if fileName != infoData_ls[i].strip(".xlsx"):
        print(
            "File Naming is wrong with either \n {0} or \n {1}".format(
                rawData_ls[i], infoData_ls[i]
            )
        )
        break
    else:
        testDate = fileName.split("_")[
            0
        ]  # extract test date, which will be stored in data_dict

        # extract participant index info
        testInfo = pd.read_excel(
            path.join(basepath, infoData_ls[i]), header=None, engine="openpyxl"
        ).values

        # load raw data and row name
        rawData = np.loadtxt(path.join(basepath, rawData_ls[i]))
        rawData_mean = np.mean(rawData, axis=1)
        # add on the following line to cap the maximum at 250
        rawData_mean = np.minimum(rawData_mean, np.ones(rawData_mean.size) * 250)

        if testInfo.size != rawData.shape[0]:
            print("Probably wrong content given in the file {}".format(infoData_ls[i]))
            break
        else:
            testInfo_split = []
            participantNo = (
                []
            )  # the list of participant no whose samples were tested in this particupar plate
            data_slice = []

            for jj in testInfo[:, 0]:
                jj = jj.strip()  # remove all possible leading and trailing whitespace
                testInfo_slice = jj.split("_")

                # make each slice has the same dimension, to convert into ndarray
                while len(testInfo_slice) < 4:
                    testInfo_slice.append(None)
                while len(testInfo_slice) > 4:
                    testInfo_slice.pop(-1)

                testInfo_split.append(testInfo_slice)

            #  import pdb;pdb.set_trace()
            # change Omicron RBD to BA1/2 for consistency and convert None to WT RBD
            for ii, Info in enumerate(testInfo_split):
                if not any('CR3022' in item if item is not None else item for item in Info):
                    if 'Omicron' in Info or 'Omicron RBD' in Info:
                        testInfo_split[ii] = \
                            ['BA1/2 RBD' if 'Omicron' in item else item for item in Info]
                    #  if  'Omicron RBD' in Info:
                        #  testInfo_split[ii] = \
                            #  ['BA1/2 RBD' if 'Omicron RBD' in item else item for item in Info]
                    elif 'WT' in Info:
                        testInfo_split[ii] = \
                            ['WT RBD' if 'WT' in item else item for item in Info]
                    elif None in Info:
                        testInfo_split[ii] = \
                            ['WT RBD' if item is None else item for item in Info]
            # change to np.array for simpler processing
            testInfo_split = np.array(testInfo_split)

            # Extract the standard signals including positive negative controls
            spike_in = []
            mask_CR3022 = ["CR3022" in item for item in testInfo_split[:, 0]]
            standard_signal = rawData_mean[mask_CR3022]
            # print(standard_signal)
            for k in testInfo_split[mask_CR3022, 0]:
                spike_in.append(
                    float(k.split(" ")[0][:-2])
                )  # extract the spike in concentrations from standards

            spike_in_EC50 = np.nan
            if (
                len(spike_in) >= 4
            ):  # consider calculating the value only when more than 3 conditions in the standard series

                # now we have [standard signal] and [the corresponding spike-in concentrations]. Calculate the spike in concentration corresponding to half maximum, i.e. spike_in_EC50
                half_max_CR3022 = (standard_signal[0] + standard_signal[-1]) / 2
                for l in range(len(spike_in) - 2):  # exclude the 0 pM CR3022 point
                    if standard_signal[l] > half_max_CR3022 > standard_signal[l + 1]:
                        y1 = standard_signal[l]
                        x1 = np.log(spike_in[l])
                        y2 = standard_signal[l + 1]
                        x2 = np.log(spike_in[l + 1])
                        spike_in_EC50 = np.exp(
                            x1 + (y1 - half_max_CR3022) / (y1 - y2) * (x2 - x1)
                        )

            # extract the participant no. from the list
            dataIndex = np.where(
                (testInfo_split[:, 1] == "Blood") | (testInfo_split[:, 1] == "Saliva")
            )  # dataIndex is an array dataIndex[0] is the list
            participantNo = list(
                set(testInfo_split[dataIndex[0], 0])
            )  # remove repetitive elements and transfer back to ordered list

            varianRBDList = list(
            set(testInfo_split[dataIndex[0], 3]))

            for m in participantNo:  # m is string
                if m.isnumeric() == True:  # just in case some element is not a number

                    # now we have str(testDate), [participantNo], array(tastInfo_split),rawData_mean. We need to fill in [data_slice], which is to be appended into [data]
                    for varianRBD in varianRBDList:
                        data_slice = np.array(
                            [np.nan] * len(data_tag)
                        )  # to hold test results for individual participant
                        index_m = np.where(
                            (testInfo_split[:, 0] == m)
                            & (testInfo_split[:, 3] == varianRBD)
                        )  # participant number matches
                        mask = np.isin(data_tag, testInfo_split[index_m, 2])
                        for n in np.where(mask == True)[0]:
                            index_m_n = np.where(
                                (testInfo_split[:, 0] == m)
                                & (testInfo_split[:, 2] == data_tag[n])
                                & (testInfo_split[:, 3] == varianRBD)
                            )

                            if len(index_m_n[0]) != 1:  # adding [0] gives list
                                print(
                                    "Repetitive test conditions found on {}".format(
                                        infoData_ls[i]
                                    )
                                )
                                raise NameError("Check the test conditions")
                            else:
                                data_slice[n] = rawData_mean[index_m_n]

                        #  breakpoint()
                        data_slice = np.insert(
                            data_slice, 0, [ m,
                                rawData_mean[
                                    np.where(testInfo_split[:, 0] == "400pM CR3022")
                                ]
                                if ("400pM CR3022" in testInfo_split[:, 0])
                                else np.amax(rawData_mean),
                                np.amin(rawData_mean),
                            ],
                        )  # rawData_mean[np.where(testInfo_split[:,0] == '0pM CR3022')] if ('0pM CR3022' in testInfo_split[:,0]) else np.amin(rawData_mean)])
                        # append, participant no, positive and negative signal in front of the test results, positive and negative will be maximum/minimum if that data is not present
                        #  import pdb; pdb.set_trace()

                        if not (
                            any(
                                np.isin(
                                    [
                                        "500pM CR3022",
                                        "400pM CR3022",
                                        "250pM CR3022",
                                        "100pM CR3022",
                                        "20pM CR3022",
                                    ],
                                    testInfo_split[:, 0],
                                )
                            )
                            and "0pM CR3022" in testInfo_split[:, 0]
                        ):
                            print("txt file {} lacks control data".format(rawData_ls[i]))

                        data_slice = np.insert(
                            data_slice, 3, (data_slice[1] + data_slice[2]) / 2
                        )  # add half maximum into the data_slice

                        data.append(data_slice)

                        #####time to calculate titer
                        blood_data = data_slice[
                            np.array(blood_sample_index) + 4
                        ]  # blood sample test results, the first 4 position has been taken by participant no., positive, negative, and half max
                        if blood_data[0] < data_slice[3] < blood_data[1]:
                            mask[
                                0
                            ] = False  # for some data the 200X is lower than half maximum and 1000X. Exclude that for accurate calculation
                        signal = (
                            [np.inf]
                            + blood_data[mask[0 : len(blood_sample_index)]].tolist()
                            + [0]
                        )  # exclude Nan from the calculation

                        dilution_factor = [
                            int(item[:-1])
                            for item in blood_sample_tag[mask[0 : len(blood_sample_index)]]
                        ]  # excluded Nan, delete the 'X' in each data label
                        for n in range(len(signal) - 1):
                            if (
                                signal[n] > data_slice[3] > signal[n + 1]
                            ):  # compare with half maximum
                                break
                        if n == 0:
                            titer = dilution_factor[0]
                        elif (
                            n == len(signal) - 2
                        ):  # if half maximum is larger than the signal from the last diution factor, use that dilution factor
                            titer = dilution_factor[-1]
                        else:

                            # now we pinpointed (the neighbour data points of the half maximum, link this two points and it crosses the half maximum. The corresponding dilution factor of the intercept is our titer)
                            # the neighbour points are
                            y1 = signal[n]
                            x1 = np.log(dilution_factor[n - 1])
                            y2 = signal[n + 1]
                            x2 = np.log(dilution_factor[n])
                            titer = np.exp(
                                x2 - (data_slice[3] - y2) / (y1 - y2) * (x2 - x1)
                            )

                        #  import pdb; pdb.set_trace()
                        # extra calculation, to arrive at the right titer with the presence of hook effect
                        if titer == 200:
                            temp = []
                            for k in range(len(signal[1:-1]) - 1):
                                if signal[k + 1] > signal[k + 2]:
                                    y1 = signal[k + 1]
                                    x1 = np.log(dilution_factor[k])
                                    y2 = signal[k + 2]
                                    x2 = np.log(dilution_factor[k + 1])
                                    intercept = np.exp(
                                        x2 - (data_slice[3] - y2) / (y1 - y2) * (x2 - x1)
                                    )
                                    temp.append(intercept)
                                else:
                                    temp.append(0)

                            if np.amax(temp) > 0:
                                titer = np.amax(temp)
                            else:
                                titer = 0

                        #####End of titer calculation

                        #####Saliva calculation
                        saliva_data = data_slice[np.array(saliva_sample_index) + 4]
                        saliva_data = np.maximum(
                            (
                                np.minimum(
                                    saliva_data,
                                    np.array([data_slice[1]] * saliva_data.size),
                                )
                                - data_slice[2]
                            )
                            / (data_slice[1] - data_slice[2]),
                            np.array([0] * saliva_data.size),
                        )  # rescale raw saliva test data according to the positive and negative signal, change negative value to zero, anything above 1 gives 1

                        titer_slice = [
                            testDate,
                            m,
                            varianRBD,
                            int(titer),
                            spike_in_EC50,
                        ] + saliva_data.tolist()
                        #####End of saliva calculation

                        #####Add all other useful information into the master file
                        usefulInfo = pa[pa["Code"] == int(m)]

                        daysInfection = -np.inf
                        if not pd.isnull(usefulInfo['Infection']).values:
                        #  if not usefulInfo['Infection'].isnull().values:
                            daysInfection = float(
                            (datetime.strptime(testDate, "%Y-%m-%d")
                            - pd.to_datetime(usefulInfo['Infection'].values[0], dayfirst=True)
                            ).days
                            )
                        #  import pdb;pdb.set_trace()
                        daysFromJabs = (
                            []
                        )  # days from the second last jab, if that exists, and days from the last jab

                        for jab in jabColNames:
                            if not pd.isnull(usefulInfo[jab]).values:
                                daysFromJabs.append(
                                    float(
                                        (
                                            datetime.strptime(testDate, "%Y-%m-%d")
                                            - pd.to_datetime(usefulInfo[jab].values[0], dayfirst=True)
                                        ).days
                                    )
                                )  # make it float so that Nan can be used later

                        while len(daysFromJabs) < len(jabColNames) + 1:
                            #  daysFromJabs.append(None)
                            #  daysFromJabs.append(np.nan)
                            daysFromJabs.append(-np.inf)
                        daysFrom2ndJab = daysFromJabs[1]
                        days1stBooster = daysFromJabs[2] # after 3rd jab
                        days2ndBooster = daysFromJabs[3]

                        mRNAList = ["Pfizer", "Moderna"]
                        vaccine = usefulInfo["Vaccine Brand"].values[0]
                        boostNum = vaccine.count("+")
                        #  vaccine = vaccine.split("+")[0]
                        vaccine = vaccine.split("+")
                        vaccine = [item.strip() for item in vaccine]
                        if all([item in mRNAList for item in vaccine]):
                            vaccineCategory = "mRNA"
                        else:
                            vaccineCategory = "mixVaccine"

                        if boostNum == 0:
                            vaccineCurrent = vaccine[0]
                        elif boostNum == 1:
                            if days1stBooster < 0:
                                vaccineCurrent = vaccine[0]
                            else:
                                vaccineCurrent = vaccine[1]
                        elif boostNum == 2:
                            if days1stBooster < 0:
                                vaccineCurrent = vaccine[0]
                            elif days1stBooster >= 0 and days2ndBooster < 0:
                                vaccineCurrent = vaccine[1]
                            elif days2ndBooster >= 0:
                                vaccineCurrent = vaccine[2]

                        #  breakpoint()


                        #  if len(daysFromJabs) < 2:
                        #      daysFromJabs = np.insert(
                        #          daysFromJabs, 0, np.nan
                            #  )  # for cases where only one jab is given

                        # for participants with multiple shots of different vaccines. The no of days is culculated based on the difference between test date and date of the first dose mRNA vaccine

                        #  elif len(daysFromJabs) == 4:
                            #  daysFromJabs.pop()

                        #  vaccine = usefulInfo["Vaccine Brand"].values[0]
                        #  if ("+" in vaccine) & (daysFromJabs[-1] <= 0):
                            #  vaccine = vaccine.split("+")[0]
                            #  daysFromJabs.pop()

                        #  daysFrom2ndJab = (
                        #      daysFromJabs[1] if "+" in vaccine else daysFromJabs[-1]
                        #  ) # after 2nd Jab

                        #  while len(daysFromJabs) < len(jabColNames) + 1:
                            #  daysFromJabs.append(None)
                        #  daysFrom2ndJab = daysFromJabs[1]
                        #  days1stBooster = daysFromJabs[2] # after 3rd jab
                        #  days2ndBooster = daysFromJabs[3]
                        titer_slice = titer_slice + [
                            vaccineCategory,
                            vaccineCurrent,
                            daysInfection,
                            daysFrom2ndJab,
                            days1stBooster,
                            days2ndBooster,
                            boostNum,
                            int(usefulInfo["Age"].values[0]),
                            usefulInfo["Gender"].values[0],
                            usefulInfo["Race"].values[0],
                            i,
                        ]
                        titerMaster.append(titer_slice)
titerMaster_df = pd.DataFrame(
    titerMaster,
    columns=[
        "Test Date",
        "Participant",
        "Varian",
        "Titer",
        "CR3022 EC50",
        "Saliva x5",
        "Saliva X20",
        "vaccineCategory",
        "vaccineCurrent",
        "daysInfection",
        "Days from 2nd Jab",
        "Days from 1st Booster",
        "Days from 2nd Booster",
        "boostNum",
        "Age",
        "Gender",
        "Race",
        "File Number",
    ],
)


pd.DataFrame(data).to_csv(
    "data.csv",
    sep="\t",
    encoding="utf-8",
    index=False,
    header=["participant", "+", "-", "half maximum"] + data_tag.tolist(),
)
#  titerMaster_df.to_csv(
    #  "titerMaster.csv",
    #  sep="\t",
    #  index=False,
    #  header=[
        #  "Test Date",
        #  "Participant",
        #  "Varian",
        #  "Titer",
        #  "CR3022 EC50",
        #  "Saliva x5",
        #  "Saliva X20",
        #  "vaccineCategory",
        #  "vaccineCurrent",
        #  "daysInfection",
        #  "Days from 2nd Jab",
        #  "Days from 1st Booster",
        #  "Days from 2nd Booster",
        #  "boostNum",
        #  "Age",
        #  "Gender",
        #  "Race",
        #  "File Number",
    #  ],
#  )
titerMaster_df.to_csv("titerMaster.csv", index=False)
print("done")
#####End of loading the useful info
