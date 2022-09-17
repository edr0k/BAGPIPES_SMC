import numpy as np
import pandas as pd
import bagpipes as pipes
from astropy.table import Table, join
import math

# Essa é a biblioteca do BAGPIPES, para baixar basta acessar: https://bagpipes.readthedocs.io/en/latest/index.html e seguir os
# os passos. O PyMultiNest é essencial e também precisa ser instalado, caso contrário o código não vai rodar, mas nesse
# mesmo link ele indica como fazer. Qualquer coisa pode me mandar email que eu te ajudo.


# A ordem dos filtros importa para o código, então tem que colocar na mesma ordem que a fotometria foi passada.
filter_list_splus = ['filters/uJAVA.dat',
                     'filters/F0378.dat',
                     'filters/F0395.dat',
                     'filters/F0410.dat',
                     'filters/F0430.dat',
                     'filters/gSDSS.dat',
                     'filters/F0515.dat',
                     'filters/rSDSS.dat',
                     'filters/F0660.dat',
                     'filters/iSDSS.dat',
                     'filters/F0861.dat',
                     'filters/zSDSS.dat']
#####################################################################################################################
# Dicionários de parâmetros fixos e livres, nessa parte existe uma boa flexibilidade do fitting.
# Segue o link da documentação do BAGPIPES, com alguns exemplos de como usar:
# https://bagpipes.readthedocs.io/en/latest/fit_instructions.html

# A primeiro momento estou usando um modelo simples de Burst, deixando a idade, metalicidade, massa de formação para
# serem ajustados. Nunca tinha usado ele e descobri que se eu faço só isso, ele não faz muito bem a SFH e não sei
# muito bem o porquê.
burst = {}
burst["age"] = (0.01, 15.)                  # Vary age from 0 to 15 Gyr
burst["metallicity"] = (0., 1.)          # Vary metallicity from 0 to 2.5 Solar (Está em Z sobre Zsol, não em [Fe/H])
burst["massformed"] = (0., 9.)           # Vary log_10(mass formed) from 0 to 9

dust = {}
dust["type"] = "Calzetti"
dust["Av"] = 0
dust["eta"] = 3

nebular = {}
nebular["logU"] = -3


fit_instructions = {}
fit_instructions["nebular"] = nebular
fit_instructions["redshift"] = (0., 0.001)
fit_instructions["dust"] = dust
fit_instructions["burst"] = burst         # Add the burst sfh component to the fit
fit_instructions["t_bc"] = 0.01           # Max age of birth clouds: Gyr

df = pd.read_csv('data/phot_100_all_derred.csv')

#Função load data necessária para o fit. Ela pega a fotometria do ID que é passado e transforma para fluxo em mJy,
# assim como a incerteza. A conta do fluxo fui eu quem fez, espero estar certa. Se não estiver, por favor me avise!
def load_data_splus(id):
        print(id)
        try:
            object, field = id.split('_')
        except:
            object1, object2, field = id.split('_')
            object = object1 + '_' + object2

        fluxes = []
        fluxerrs = []
        galaxy_param = df[(df['object'] == object) & (df['field'] == field)]

        for k in range(5, 5+len(filter_list_splus)):
            m = galaxy_param[galaxy_param.columns[k]].values
            if math.isnan(m) or m == np.inf or m == 99.:
                f = np.array([99.])
                delta_f = np.array([99.])
            else:
                f = 10**(9.56) * 10**(-m/2.5)  # flux in mJy
                delta_m = galaxy_param[galaxy_param.columns[12 + k]].values
                delta_f = f * (1/2.5) * np.log(10) * delta_m

            #print(m, delta_m, f, delta_f)

            fluxes.append(f)  #mJy
            fluxerrs.append(delta_f)  #mJy
        photometry = np.c_[fluxes, fluxerrs]

        print(photometry)
        return photometry

IDs = df['object'] + '_'+ df['field']
#redshifts = df['z']
fit_cat = pipes.fit_catalogue(IDs, fit_instructions, load_data_splus, spectrum_exists=False,
                              cat_filt_list=filter_list_splus, run="phot_100_all_derred_burst", make_plots=True,
                              full_catalogue=True)

fit_cat.fit(verbose=False)

results_cat = Table.read('pipes/cats/phot_100_all_derred_burst.fits', format='fits')
literature_cat = Table.read('literature/parameter-cat_Conly_ref.csv', format='csv')

from astropy.table import hstack
results_cat = hstack([results_cat, literature_cat])
results_cat.write('pipes/cats/phot_100_all_derred_burst.csv', format='csv', overwrite=True)

