import cobra
import numpy as np
from modules.zavrel_FIT import photodamage_helper, fitted_glyc, proteincontent
from modules.add_glycogen import update_prot_glyc

blue_light = 27.5


class Sensitivity_Calculator:
    def __init__(
        self,
        model: cobra.Model,
        kd: float,
        kl: float,
        alpha: float,
        alpha2: float,
        maint: float
    ):
        self.model = model
        self.kd = kd
        self.kl = kl
        self.alpha = alpha
        self.alpha2 = alpha2
        self.maint = maint

        self.params = np.array([
            self.kl,
            self.kd,
            self.alpha,
            self.alpha2
        ])
        self.default_mu = {}

    def calc_sens(self, point: float, verbose: bool = False) -> dict:

        # calculate default growth rates and propagate dict:
        with self.model:
            self.defmus = photodamage_helper(
                self.model,
                'EX_E1_ext_b',
                "BM0009",
                -(point + blue_light),
                self.kl,
                self.kd,
                self.alpha,
                self.alpha2,
                False,
                (True, 2),
                self.maint
            )
        self.defmus

        self.sensitivities = {
            "kl": '',
            "kd": '',
            "a": '',
            "a2": '',
            "I0": ''
        }
        self.params_lb = self.params * .99
        self.params_ub = self.params * 1.01

        # make list of lists with one parameter adjusted per list:
        self.ub_params_list = [
            [
                self.params[ind] if ind != i
                else self.params_ub[i]
                for ind in range(len(self.params))
            ] + [point]
            for i in range(len(self.params))
        ] + [[i for i in self.params] + [(point) * 1.01]]
        self.lb_params_list = [
            [
                self.params[ind] if ind != i
                else self.params_lb[i]
                for ind in range(len(self.params))
            ] + [point]
            for i in range(len(self.params))
        ] + [[i for i in self.params] + [(point) * .99]]

        # calculate growth rates for every set of parameters (but light)
        for i in range(len(self.sensitivities)):
            with self.model:
                self.sensitivities[list(self.sensitivities)[i]] = (
                    (
                        (mu_ub := photodamage_helper(
                            update_prot_glyc(
                                self.model,
                                proteincontent(point + blue_light),
                                fitted_glyc(point + blue_light)
                            ),
                            'EX_E1_ext_b',
                            "BM0009",
                            -(self.ub_params_list[i][-1] + blue_light),
                            self.ub_params_list[i][0],
                            self.ub_params_list[i][1],
                            self.ub_params_list[i][2],
                            self.ub_params_list[i][3],
                            0,
                            (False, 2),
                            self.maint
                        ))
                        -
                        (mu_lb := photodamage_helper(
                            update_prot_glyc(
                                self.model,
                                proteincontent(point + blue_light),
                                fitted_glyc(point + blue_light)
                            ),
                            'EX_E1_ext_b',
                            "BM0009",
                            -(self.lb_params_list[i][-1] + blue_light),
                            self.lb_params_list[i][0],
                            self.lb_params_list[i][1],
                            self.lb_params_list[i][2],
                            self.lb_params_list[i][3],
                            0,
                            (False, 2),
                            self.maint
                        ))
                    )
                    /
                    self.defmus
                    /
                    (
                        (self.ub_params_list[i][i] - self.lb_params_list[i][i])
                        /
                        ([j for j in self.params] + [point])[i]
                    )
                )
                if verbose:
                    print(self.defmus)
                    print(
                        list(self.sensitivities)[i],
                        (list(self.params)+[point])[i],
                        self.lb_params_list[i][i],
                        self.ub_params_list[i][i],
                        mu_lb, mu_ub
                    )
        return self.sensitivities
