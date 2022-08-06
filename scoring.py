import freq

dfAFR = freq.dfAFR
dfAMR = freq.dfAMR
dfEAS = freq.dfEAS
dfEUR = freq.dfEUR
dfSAS = freq.dfSAS

power = {'rs1815739': {'T|T':0,'T|C':1,'C|C':2}, 'rs4253778': {'G|G':0, 'C|G':1, 'C|C':2}, 'rs1042713': {'A|A':0, 'G|A':1, 'G|G':2},
         'rs8192678': {'C|C':0, 'T|C':1, 'T|T':2}, 'rs2070744': {'C|C':0, 'T|C':1, 'T|T':2}, 'rs11549465': {'C|C':0, 'T|C': 1, 'T|T':2}}

endurance = {'rs1815739': {'T|T':2,'T|C':1,'C|C':0}, 'rs4253778': {'G|G':2, 'C|G':1, 'C|C':0}, 'rs1042713': {'A|A':2, 'G|A':1, 'G|G':0},
             'rs8192678': {'C|C':2, 'T|C':1, 'T|T':0}, 'rs2070744': {'C|C':2, 'T|C':1, 'T|T':0}, 'rs11549465': {'C|C':0, 'T|C': 1, 'T|T':2}}

injuryrisk = {'rs1800012': {'C|C':2, 'C|A':1, 'A|A':0}, 'rs12722': {'T|T':2,'T|C':1,'C|C':0}, 'rs1049434': {'A|A':2, 'T|A':1, 'T|T':0},
              'rs1800795': {'C|C':2, 'C|G':1, 'G|G':0}}

oxygencapacity = {'rs2070744': {'C|C':2, 'T|C':1, 'T|T':0}, 'rs11549465': {'C|C':2, 'T|C': 1, 'T|T':0}}

motorskills = { 'rs6265': {'C|C':2, 'T|C':1, 'T|T':0}, 'rs4680': {'A|A':2, 'G|A':1, 'G|G':0}}

snpname = {'rs1815739':0, 'rs4253778':2, 'rs1042713':3, 'rs8192678':4, 'rs2070744':5, 'rs11549465':6, 'rs1800012':7, 'rs12722':8,
           'rs1049434':9, 'rs1800795':10, 'rs6265':11, 'rs4680':12}

def geneticscore(df):
    powerscr = []
    endurancescr = []
    injuryriskscr = []
    oxygencapacityscr = []
    motorskillsscr = []

    forpower = [0, 2, 3, 4, 5, 6]
    forendurance = [0, 2, 3, 4, 5, 6]
    forinjuryrisk = [7, 8, 9, 10]
    foroxygencapacity = [5, 6]
    formotorskills = [11, 12]

    dfpower = df.iloc[:, forpower]
    dfpower = dfpower.replace(power)
    powerscr.append((dfpower.sum(axis=1) / (2 * len(forpower))) * 100)

    dfendurance = df.iloc[:, forendurance]
    dfendurance = dfendurance.replace(endurance)
    endurancescr.append((dfendurance.sum(axis=1) / (2 * len(forendurance))) * 100)

    dfinjuryrisk = df.iloc[:, forinjuryrisk]
    dfinjuryrisk = dfinjuryrisk.replace(injuryrisk)
    injuryriskscr.append((dfinjuryrisk.sum(axis=1) / (2 * len(forinjuryrisk))) * 100)

    dfoxygencapacity = df.iloc[:, foroxygencapacity]
    dfoxygencapacity = dfoxygencapacity.replace(oxygencapacity)
    oxygencapacityscr.append((dfoxygencapacity.sum(axis=1) / (2 * len(foroxygencapacity))) * 100)

    dfmotorskills = df.iloc[:, formotorskills]
    dfmotorskills = dfmotorskills.replace(motorskills)
    motorskillsscr.append((dfmotorskills.sum(axis=1) / (2 * len(formotorskills))) * 100)

    return powerscr, endurancescr, injuryriskscr, oxygencapacityscr, motorskillsscr

powerscr_AFR, endurancescrAFR, injuryriskscrAFR, oxygencapacityscrAFR, motorskillsscrAFR = geneticscore(dfAFR)
powerscr_AMR, endurancescrAMR, injuryriskscrAMR, oxygencapacityscrAMR, motorskillsscrAMR = geneticscore(dfAMR)
powerscr_EAS, endurancescrEAS, injuryriskscrEAS, oxygencapacityscrEAS, motorskillsscrEAS = geneticscore(dfEAS)
powerscr_EUR, endurancescrEUR, injuryriskscrEUR, oxygencapacityscrEUR, motorskillsscrEUR = geneticscore(dfEUR)
powerscr_SAS, endurancescrSAS, injuryriskscrSAS, oxygencapacityscrSAS, motorskillsscrSAS = geneticscore(dfSAS)