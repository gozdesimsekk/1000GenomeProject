import pandas as pd
import matplotlib.pyplot as plt
import numpy as np
from scipy.stats import chisquare
from tabulate import tabulate

dfAFR = pd.concat(map(pd.read_csv, ['output/rs1815739/AFR.csv', 'output/rs4644994/AFR.csv', 'output/rs4253778/AFR.csv', 'output/rs1042713/AFR.csv', 'output/rs8192678/AFR.csv', 'output/rs2070744/AFR.csv',
                                 'output/rs11549465/AFR.csv','output/rs1800012/AFR.csv', 'output/rs12722/AFR.csv', 'output/rs1049434/AFR.csv', 'output/rs1800795/AFR.csv', 'output/rs6265/AFR.csv', 'output/rs4680/AFR.csv']), ignore_index=True, axis=1)
dfAFR = dfAFR.iloc[:,[1, 6 , 11, 16, 21, 26, 31, 36, 41, 46, 51, 56, 61]]

dfAMR = pd.concat(map(pd.read_csv, ['output/rs1815739/AMR.csv', 'output/rs4644994/AMR.csv', 'output/rs4253778/AMR.csv', 'output/rs1042713/AMR.csv', 'output/rs8192678/AMR.csv', 'output/rs2070744/AMR.csv',
                                 'output/rs11549465/AMR.csv','output/rs1800012/AMR.csv', 'output/rs12722/AMR.csv', 'output/rs1049434/AMR.csv' ,'output/rs1800795/AMR.csv', 'output/rs6265/AMR.csv', 'output/rs4680/AMR.csv']), ignore_index=True, axis=1)
dfAMR = dfAMR.iloc[:,[1, 6 , 11, 16, 21, 26, 31, 36, 41, 46, 51, 56, 61]]

dfEAS = pd.concat(map(pd.read_csv, ['output/rs1815739/EAS.csv', 'output/rs4644994/EAS.csv', 'output/rs4253778/EAS.csv', 'output/rs1042713/EAS.csv', 'output/rs8192678/EAS.csv', 'output/rs2070744/EAS.csv',
                                 'output/rs11549465/EAS.csv','output/rs1800012/EAS.csv', 'output/rs12722/EAS.csv', 'output/rs1049434/EAS.csv' ,'output/rs1800795/EAS.csv', 'output/rs6265/EAS.csv', 'output/rs4680/EAS.csv']), ignore_index=True, axis=1)
dfEAS = dfEAS.iloc[:,[1, 6 , 11, 16, 21, 26, 31, 36, 41, 46, 51, 56, 61]]

dfEUR = pd.concat(map(pd.read_csv, ['output/rs1815739/EUR.csv', 'output/rs4644994/EUR.csv', 'output/rs4253778/EUR.csv', 'output/rs1042713/EUR.csv', 'output/rs8192678/EUR.csv', 'output/rs2070744/EUR.csv',
                                 'output/rs11549465/EUR.csv','output/rs1800012/EUR.csv', 'output/rs12722/EUR.csv', 'output/rs1049434/EUR.csv' ,'output/rs1800795/EUR.csv', 'output/rs6265/EUR.csv', 'output/rs4680/EUR.csv']), ignore_index=True, axis=1)
dfEUR = dfEUR.iloc[:,[1, 6 , 11, 16, 21, 26, 31, 36, 41, 46, 51, 56, 61]]

dfSAS = pd.concat(map(pd.read_csv, ['output/rs1815739/SAS.csv', 'output/rs4644994/SAS.csv', 'output/rs4253778/SAS.csv', 'output/rs1042713/SAS.csv', 'output/rs8192678/SAS.csv', 'output/rs2070744/SAS.csv',
                                 'output/rs11549465/SAS.csv','output/rs1800012/SAS.csv', 'output/rs12722/SAS.csv', 'output/rs1049434/SAS.csv' ,'output/rs1800795/SAS.csv', 'output/rs6265/SAS.csv', 'output/rs4680/SAS.csv']), ignore_index=True, axis=1)
dfSAS = dfSAS.iloc[:,[1, 6 , 11, 16, 21, 26, 31, 36, 41, 46, 51, 56, 61]]

for i in [dfAFR, dfAMR, dfEAS, dfEUR, dfSAS]:
    i.columns = ["rs1815739", "rs4644994", "rs4253778", "rs1042713", "rs8192678", "rs2070744", "rs11549465", "rs1800012", "rs12722", "rs1049434","rs1800795", "rs6265", "rs4680"]
cols = ["rs1815739", "rs4644994", "rs4253778", "rs1042713", "rs8192678", "rs2070744", "rs11549465", "rs1800012", "rs12722",
     "rs1049434", "rs1800795", "rs6265", "rs4680"]
#etnik kökenlerde  heterozigotları tek fenotipte göstermek için değişim
def converter(df):
    heterovals = [["A|C","C|A"],["A|T","T|A"],["A|G","G|A"],["C|T","T|C"],["G|C","C|G"],["T|G","G|T"]]
    for i in range(len(df)):
        for j in range(len(cols)):
            for k,l in heterovals:
                if k == df.iloc[i,j]:
                    df.iloc[i,j] = l
    return df
dfAFR, dfAMR, dfEAS, dfEUR, dfSAS = converter(dfAFR), converter(dfAMR), converter(dfEAS), converter(dfEUR), converter(dfSAS)
populations = ["AFR", "AMR", "EAS", "EUR", "SAS"]
#genotype frekanslarının hesaplanması
gozlemlenen = []
beklenen = []
dfERRORS = {}
for index, i in enumerate([dfAFR, dfAMR, dfEAS, dfEUR, dfSAS]):
    for j in cols:
        count = i[j].value_counts()
        dfcount = pd.DataFrame(count)
        N = dfcount.sum() #toplam genotip sayısı
        freq = (dfcount / N)
#allele frekanslarının hesaplanması
        list = i.filter([j]).values.reshape(1, -1).ravel().tolist()
        unsymbol = []
        symbol = '|'
        for element in list:
            temp = ""
            for ch in element:
                if ch not in symbol:
                    temp += ch
            unsymbol.append(temp)
        splitted = [k for b in unsymbol for k in b]
        count = dict(map(lambda x: (x, splitted.count(x)), splitted))
        freqofallels = {key: val / sum(count.values()) for key,
                                                                 val in count.items()}
#beklenen genotip frekanslarının hesaplanması
        minorallele = (min(freqofallels.values()))
        minoralleletype = min(count.items(), key=lambda x: x[1])
        h0genotype = ((1-minorallele) ** 2)
        h1genotype = 2 * (minorallele * (1-minorallele))
        h2genotype = (minorallele ** 2)
        freq = freq.to_dict()
        dict(sorted(freq.items(), key=lambda item: item[1]))
#beklenen genotip frekanslarının listelenmesi(çoktan aza)
        estimatedgenotypes = [h0genotype, h1genotype, h2genotype]
        estimatedgenotypes.sort(reverse=True)
#görülmeyen gözlemlenen genotip frekanslarına 3. frekans olarak 0 ekleme.
        listofgenotypes = []
        listoffreq = []
        for k, v in freq.items():
            for k1, v1 in v.items():
                listoffreq.append(v1)
                listofgenotypes.append(k1)
            if len(listoffreq) == 2:
                listoffreq.append(0)
        genotypes = []
        symbol = '|'
        for element in listofgenotypes:
            temp = ""
            for ch in element:
                if ch not in symbol:
                    temp += ch
            genotypes.append(temp)
#görülmeyen genotipler için görülmeyen genotipi ekleme
        if len(genotypes) < 3:
            for gent in genotypes:
                if gent[0] != gent[1]:
                    if "{}{}".format(gent[0],gent[0]) not in genotypes:
                        genotypes.append("{}{}".format(gent[0],gent[0]))
                    elif "{}{}".format(gent[1],gent[1]) not in genotypes:
                        genotypes.append("{}{}".format(gent[1], gent[1]))
#hata hesaplaması ve görselleştirme
        error = np.divide(np.square(np.subtract(listoffreq,estimatedgenotypes)), estimatedgenotypes)
        errorsum = error.sum()
        errorl = error.tolist()
        dferror = pd.DataFrame(errorl,genotypes, columns=[j])
        ax = dferror.T.plot(kind='barh', color=['pink', 'cyan', 'crimson'])
        plt.xlabel(" error ")
        plt.title(populations[index])
        plt.gcf().set_size_inches(10, 5)
        plt.close()
        dr = dferror.to_dict()
        dfERRORS.update(dr)
        #güven aralığı hesaplama %95
        arrayy = np.divide(errorsum, (np.array(N) ** 0.5))
        up = estimatedgenotypes + (1.96 * arrayy)
        low = estimatedgenotypes - (1.96 * arrayy)
        confidence_interval = []
        for s in range(len(up)):
            confidence_interval.append([low[s], up[s]])
        CI = dict(zip(genotypes, confidence_interval))

#kikare testi
        #obs ve exp birey sayısını hesaplamak için
        for q in listoffreq:
            gozlemlenen.append(q * N.values[0])
        for w in estimatedgenotypes:
            beklenen.append(w * N.values[0])
#hataları excelde göstermek için

'''dfERRORS = pd.DataFrame(dfERRORS)
dfERRORS.to_excel('SNPerror.xlsx')'''

#etnik kökenlere göre countları ayırmak
obslist = np.array_split(gozlemlenen, len(gozlemlenen)/3)
AFRobs = obslist[0:13]
AMRobs = obslist[13:26]
EASobs = obslist[26:39]
EURobs = obslist[39:52]
SASobs = obslist[52:65]

obsAFRdict = dict(zip(cols, AFRobs))
obsAMRdict = dict(zip(cols, AMRobs))
obsEASdict = dict(zip(cols, EASobs))
obsEURdict = dict(zip(cols, EURobs))
obsSASdict = dict(zip(cols, SASobs))
def adict(dict1,dict2):
    for i in dict1.keys():
        for j in dict2.keys():
            if i == j:
                dict1[i] = dict1[i] + dict2[j]
    for j in dict2.keys():
        if j not in dict1.keys():
            a = {j : dict2[j]}
            dict1.update(a)
    return dict1
tot1 = adict(obsAFRdict, obsAMRdict)
tot2 = adict(tot1, obsEASdict)
tot3 = adict(tot2, obsEURdict)
observed = adict(tot3, obsSASdict)

explist = np.array_split(beklenen, len(beklenen)/3)
AFRexp = explist[0:13]
AMRexp = explist[13:26]
EASexp = explist[26:39]
EURexp = explist[39:52]
SASexp = explist[52:65]
expAFRdict = dict(zip(cols, AFRexp))
expAMRdict = dict(zip(cols, AMRexp))
expEASdict = dict(zip(cols, EASexp))
expEURdict = dict(zip(cols, EURexp))
expSASdict = dict(zip(cols, SASexp))
def adict(dict1,dict2):
    for i in dict1.keys():
        for j in dict2.keys():
            if i == j:
                dict1[i] = dict1[i] + dict2[j]
    for j in dict2.keys():
        if j not in dict1.keys():
            a = {j : dict2[j]}
            dict1.update(a)
    return dict1
tott1 = adict(expAFRdict, expAMRdict)
tott2 = adict(tott1, expEASdict)
tott3 = adict(tott2, expEURdict)
expected = adict(tott3, expSASdict)

listt = []
for l in cols:
    listt.append(chisquare(observed[l], f_exp=expected[l], ddof=1))
dcp = dict(zip(cols,listt))
idx = ['chi score', 'p value']
pvaltable = tabulate(dcp, headers=cols, tablefmt="fancy_grid", showindex=idx)
#print(pvaltable)
headers = ["rs_id", "errors"]
errortable = tabulate(dfERRORS.items(), headers = headers)
#print(errortable)










