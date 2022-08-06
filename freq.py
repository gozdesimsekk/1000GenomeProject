import pandas as pd
import matplotlib.pyplot as plt

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

def converter(df):
    heterovals = [["A|C","C|A"],["A|T","T|A"],["A|G","G|A"],["C|T","T|C"],["G|C","C|G"],["T|G","G|T"]]
    for i in range(len(df)):
        for j in range(len(cols)):
            for k,l in heterovals:
                if k == df.iloc[i,j]:
                    df.iloc[i,j] = l
    return df
dfAFR, dfAMR, dfEAS, dfEUR, dfSAS = converter(dfAFR), converter(dfAMR), converter(dfEAS), converter(dfEUR), converter(dfSAS)
populations = ["AFR", "AMR","EAS", "EUR", "SAS"]
for index, i in enumerate([dfAFR, dfAMR, dfEAS, dfEUR, dfSAS]):
    for j in cols:
        count = i[j].value_counts()
        dfcount = pd.DataFrame(count)
        total = dfcount.sum()
        freq = (dfcount / total) * 100
        print(freq)
        dffreq = pd.DataFrame(freq)
        ax = dffreq.T.plot(kind='barh', stacked=True, color = ['pink', 'cyan', 'crimson'])
        for p in ax.patches:
            h, w, x, y = p.get_height(), p.get_width(), p.get_x(), p.get_y()
            text = f'{w :0.2f} %'
            ax.annotate(text=text, xy=(x + w / 2, y + h / 2), ha='left', va='center_baseline', color='black', size=15)
        plt.xlabel(" % Frequency")
        plt.title(populations[index])
        plt.gcf().set_size_inches(20, 5)
        plt.savefig("output/genotypes/{}_{}.png".format(populations[index], j), dpi=200)
        plt.close()

populations = ["AFR", "AMR", "EAS", "EUR", "SAS"]

for index, i in enumerate([dfAFR, dfAMR, dfEAS, dfEUR, dfSAS]):
    for j in cols:
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
        freqofallels = {key: val / sum(count.values()) *100 for key,
                                                  val in count.items()}
        
        allel_df = pd.DataFrame(freqofallels, index=[j])
        af = allel_df.T.plot(kind='barh', color='darkviolet')
        for p in af.patches:
            h, w, x, y = p.get_height(), p.get_width(), p.get_x(), p.get_y()
            text = f'{w :0.2f} %'
            af.annotate(text=text, xy=(x + w / 2, y + h / 2), ha='left', va='center_baseline', color='black', size=15)
        plt.xlabel("% Frequency")
        plt.ylabel("Allele Types")
        plt.title(populations[index])
        plt.gcf().set_size_inches(20, 5)
        plt.savefig("output/allels/{}_{}.png".format(populations[index], j), dpi=200)
        plt.close()


