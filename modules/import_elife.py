import pandas as pd
import numpy as np
import uncertainties
import matplotlib.pyplot as plt


class DataImporter:

    def __init__(self, path="zavrel/elife-42508-fig2-data1-v2.xlsx"):
        self.path = path
        self.ALL = pd.read_excel(self.path, usecols="C:J", nrows=236)
        self.DATA = {
            "light": self.ALL.iloc[0].T.to_numpy(),
            "growth": self.ALL.iloc[12].T.to_numpy(),
            "growth error": self.ALL.iloc[13].T.to_numpy(),
            "DW": self.ALL.iloc[25].T.to_numpy(),
            "DW error": self.ALL.iloc[26].T.to_numpy(),
            "light cap": self.ALL.iloc[64].T.to_numpy(),
            "light cap error": self.ALL.iloc[65].T.to_numpy(),
            "O2": self.ALL.iloc[77].T.to_numpy(),
            "O2 error": self.ALL.iloc[78].T.to_numpy(),
            "CO2": self.ALL.iloc[90].T.to_numpy(),
            "CO2 error": self.ALL.iloc[91].T.to_numpy(),
            "dark resp": self.ALL.iloc[103].T.to_numpy(),
            "dark resp error": self.ALL.iloc[104].T.to_numpy(),
            "chl": self.ALL.iloc[142].T.to_numpy(),
            "chl error": self.ALL.iloc[143].T.to_numpy(),
            "glycogen": self.ALL.iloc[194].T.to_numpy(),
            "glycogen error": self.ALL.iloc[195].T.to_numpy(),
            "protein": self.ALL.iloc[207].T.to_numpy(),
            "protein error": self.ALL.iloc[208].T.to_numpy(),
            "carbon": self.ALL.iloc[116].T.to_numpy(),
            "carbon error": self.ALL.iloc[117].T.to_numpy(),
            "nitrogen": self.ALL.iloc[129].T.to_numpy(),
            "nitrogen error": self.ALL.iloc[130].T.to_numpy()
        }

        self.DATA["light error"] = np.zeros(len(self.DATA["light"]))

        for key in self.DATA.copy():
            # add value, error pairs:
            try:
                self.DATA[key + " with error"] = np.array(
                    [
                        uncertainties.ufloat(gr, er)
                        for gr, er
                        in zip(self.DATA[key], self.DATA[key + " error"])
                    ]
                )
            except KeyError:
                continue

        self.DATA["light calc"] = (
            self.DATA["light cap"] / (self.DATA["DW"] * 0.02) * (60*60) / 1000
        )
        self.DATA["light calc with error"] = (
            self.DATA["light cap with error"] /
            (self.DATA["DW with error"] * 0.02) *
            (60*60) / 1000
        )

        self.DATA["net O2 with error"] = (
            self.DATA["O2 with error"] + self.DATA["dark resp with error"]
        )  # * 1.0324
        self.DATA["O2 fba with error"] = self.conv_for_fba(
            self.DATA["O2 with error"]
        )
        self.DATA["dark resp fba with error"] = self.conv_for_fba(
            self.DATA["dark resp with error"]
        )
        self.DATA["CO2 fba with error"] = self.conv_for_fba(
            self.DATA["CO2 with error"]
        )
        self.DATA["net O2 fba with error"] = self.conv_for_fba(
            self.DATA["net O2 with error"]
        )

        self.ERROR = {
            key[:-11]: self.DATA[key]
            for key in self.DATA
            if " with error" in key
        }

    def conv_for_fba(self, input: np.array) -> np.array:
        '''
        Convert units from per g chlorophyll to per g dry weight.

        This function converts units of a given input array from a measurement
        per gram of chlorophyll to a measurement per gram
        of dry weight. It uses known values of chlorophyll
        content and dry weight to perform the conversion.

        Parameters:
        input (np.array): An array containing values in per g chlorophyll.

        Returns:
        np.array: An array with values converted to per g dry weight.

        Notes:
        - The conversion factors are based on the provided chlorophyll content
        (CHL) and dry weight (DW) values.
        - The input values are assumed to be in micromoles (µmol) of gas per
        gram of chlorophyll per hour.
        - The output values are in millimoles (mmol) of gas per gram of dry
        weight per hour.
        '''
        M_CHL = 893.51  # g/mol chlorophyll
        CHL = self.DATA["chl with error"]
        DW = self.DATA["DW with error"]

        mgChl_per_mgDW = CHL/DW
        µmolgas_per_mgChl = input/(M_CHL*1000)
        µmolgas_per_mgDW = µmolgas_per_mgChl * mgChl_per_mgDW
        mmolgas_per_gDW_per_hour = µmolgas_per_mgDW * (60*60) * 1000

        return mmolgas_per_gDW_per_hour

    def get(self, key: str) -> np.array:
        return self.DATA[key]

    def value(self, key: str) -> np.array:
        nominals = np.array(
            [i.n for i in self.ERROR[key]]
        )
        return nominals

    def error(self, key: str) -> np.array:
        errors = np.array(
            [i.s for i in self.ERROR[key]]
        )
        return errors


class Import_Theune():
    def __init__(self):
        self.X = np.array([
            63.29787234042331,
            112.23404255318749,
            168.61702127658978,
            230.319148936162,
            316.4893617021164
        ])
        self.Y = np.array([
            0.3417721518987342,
            0.3063291139240507,
            0.3265822784810126,
            0.3493670886075948,
            0.3721518987341772
        ])
        # reproduce the calculation for ATP/NADPH as in Theune et al.:
        self.ratio = (self.Y/(1-self.Y) * 4 * 4 / 4.66 + 12 / 4.66) / 2


if __name__ == "__main__":
    zavrel = DataImporter()
    print(zavrel.get("DW with error"))
    print(zavrel.get("light cap with error"))
    print(zavrel.get("light calc"))
    print(zavrel.get("light calc with error"))
    print(zavrel.error("light calc"))
    print(f"At {zavrel.get('light')[5]} nm (gr: {zavrel.get('growth')[5]}):")
    print("Protein content:", zavrel.get("protein")[5]/zavrel.get("DW")[5])
    print("Glycogen content:", zavrel.get("glycogen")[5]/zavrel.get("DW")[5])

    fig, ax = plt.subplots()
    ax.errorbar(
        zavrel.get("light calc"), zavrel.value("O2 fba"),
        fmt="ko",
        xerr=zavrel.error("light calc"),
        yerr=zavrel.error("O2 fba"),
        capsize=2,
        elinewidth=.5
    )
    # fig.savefig("/home/hoeper/tempfiles/errorplot.svg")
    theune = Import_Theune()
    print(theune.ratio)
