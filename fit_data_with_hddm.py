import os
import hddm
from patsy import dmatrix  # for generation of (regression) design matrices
import numpy as np         # for basic matrix operations
import matplotlib.pyplot as plt
import seaborn as sns
import kabuki

model_dir = os.getcwd()

### Load Data ###
#everything should be comma seperated
data = hddm.load_csv("hddmAllData.csv")

samples = 100000
data['response'] = data['cor'] #accuracy coding
data['coh'] = data['coh'].astype(object) #convert to factor

#############################################
#Fit the model in which both drift and bound depend on coherence and volatility
m = hddm.HDDM(data, depends_on={'v': {'coh', 'vc'}, 'a': {'coh','vc'}}, p_outlier=.05)
m.find_starting_values()
m.sample(samples, burn=samples/10, thin=2, dbname=os.path.join(model_dir, 'fit_hddm_traces_1'), db='pickle')
m.save(os.path.join(model_dir, 'fit_hddm_1'))

run = False
if run:

    ##############################################
    #Check the output, note this code assumes you ran the model three times, naming it fit_hddm_1, fit_hddm_2 and fit_hddm_3
    models=[]
    for i in range(1,4):
        models.append(hddm.load('fit_hddm_%s' %i))
    from kabuki.analyze import gelman_rubin
    gelman_rubin(models)
    m = kabuki.utils.concat_models(models)

    # point estimates:
    results = m.gen_stats()
    results.to_csv(os.path.join(model_dir, 'hddm_estimates.csv'))

    # analytic plots:
    m.plot_posterior_predictive() #blue:posterior predictive;red=rt
    m.plot_posteriors()

    #Simulate data from our model and check convergence with empirical data
    ppc_data = hddm.utils.post_pred_gen(m)
    ppc_compare = hddm.utils.post_pred_stats(data, ppc_data)
    print ppc_compare

    #Grabv posteriors of bound, drift and sigma (==1)
    a00, a01, a050, a051, a10, a11, a20, a21, a40, a41 = m.nodes_db.node[['a(0.0.0.0)','a(0.0.1.0)','a(0.05.0.0)','a(0.05.1.0)','a(0.1.0.0)','a(0.1.1.0)','a(0.2.0.0)','a(0.2.1.0)','a(0.4.0.0)','a(0.4.1.0)']]
    v00, v01, v050, v051, v10, v11, v20, v21, v40, v41 = m.nodes_db.node[['v(0.0.0.0)','v(0.0.1.0)','v(0.05.0.0)','v(0.05.1.0)','v(0.1.0.0)','v(0.1.1.0)','v(0.2.0.0)','v(0.2.1.0)','v(0.4.0.0)','v(0.4.1.0)']]
    s00 = np.ones(135000);s01=s00;s050=s00;s051=s00;s10=s00;s11=s00;s20=s00;s21=s00;s40=s00;s41=s00;

    #Scale these until the bound is fixed to 1
    s00 = s00/a00.trace();v00 = v00.trace()/a00.trace();a00 = a00.trace()/a00.trace()
    s01 = s01/a01.trace();v01 = v01.trace()/a01.trace();a01 = a01.trace()/a01.trace()
    s050 = s050/a050.trace();v050 = v050.trace()/a050.trace();a050 = a050.trace()/a050.trace()
    s051 = s051/a051.trace();v051 = v051.trace()/a051.trace();a051 = a051.trace()/a051.trace()
    s10 = s10/a10.trace();v10 = v10.trace()/a10.trace();a10 = a10.trace()/a10.trace()
    s11 = s11/a11.trace();v11 = v11.trace()/a11.trace();a11 = a11.trace()/a11.trace()
    s20 = s20/a20.trace();v20 = v20.trace()/a20.trace();a20 = a20.trace()/a20.trace()
    s21 = s21/a21.trace();v21 = v21.trace()/a21.trace();a21 = a21.trace()/a21.trace()
    s40 = s40/a40.trace();v40 = v40.trace()/a40.trace();a40 = a40.trace()/a40.trace()
    s41 = s41/a41.trace();v41 = v41.trace()/a41.trace();a41 = a41.trace()/a41.trace()

    #CHECK DRIFTS
    sns.distplot(v00, hist = False);sns.distplot(v01, hist = False)
    sns.distplot(v050, hist = False);sns.distplot(v051, hist = False)
    sns.distplot(v10, hist = False);sns.distplot(v11, hist = False)
    sns.distplot(v20, hist = False);sns.distplot(v21, hist = False)
    sns.distplot(v40, hist = False);sns.distplot(v41, hist = False)
    #main effects of coherence (independent of variance)
    (((v00+v01)/2)>((v050+v051)/2)).mean()
    (((v050+v051)/2)>((v10+v11)/2)).mean()
    (((v10+v11)/2)>((v20+v21)/2)).mean()
    (((v20+v21)/2)>((v40+v41)/2)).mean()
    #no effect of variance, per coherence level
    (v00 < v01).mean()
    (v050 > v051).mean()
    (v10 > v11).mean()
    (v20 > v21).mean()
    (v40 < v41).mean()

    #CHECK SIGMA
    sns.distplot(s00, hist = False);sns.distplot(s01, hist = False)
    sns.distplot(s050, hist = False);sns.distplot(s051, hist = False)
    sns.distplot(s10, hist = False);sns.distplot(s11, hist = False)
    sns.distplot(s20, hist = False);sns.distplot(s21, hist = False)
    sns.distplot(s40, hist = False);sns.distplot(s41, hist = False)
    #No main effects of coherence (independent of variance)
    (((s00+s01)/2)>((s050+s051)/2)).mean()
    (((s050+s051)/2)<((s10+s11)/2)).mean()
    (((s10+s11)/2)<((s20+s21)/2)).mean()
    (((s20+s21)/2)>((s40+s41)/2)).mean()
    #Main effect of variance, independent from coherence
    sns.distplot(((s00+s050+s10+s20+s40)/5),hist=False)
    sns.distplot(((s01+s051+s11+s21+s41)/5),hist=False)
    (((s00+s050+s10+s20+s40)/5)>((s01+s051+s11+s21+s41)/5)).mean() #p = .014
    #Main effect of variance, per coherence level
    (s00 > s01).mean()
    (s050 > s051).mean()
    (s10 > s11).mean()
    (s20 > s21).mean()
    (s40 > s41).mean()

    #save these traces and plot in R
    np.savetxt("v00.csv", v00, delimiter=",")
    np.savetxt("v01.csv", v01, delimiter=",")
    np.savetxt("v050.csv", v050, delimiter=",")
    np.savetxt("v051.csv", v051, delimiter=",")
    np.savetxt("v10.csv", v10, delimiter=",")
    np.savetxt("v11.csv", v11, delimiter=",")
    np.savetxt("v20.csv", v20, delimiter=",")
    np.savetxt("v21.csv", v21, delimiter=",")
    np.savetxt("v40.csv", v40, delimiter=",")
    np.savetxt("v41.csv", v41, delimiter=",")

    np.savetxt("s00.csv", s00, delimiter=",")
    np.savetxt("s01.csv", s01, delimiter=",")
    np.savetxt("s050.csv", s050, delimiter=",")
    np.savetxt("s051.csv", s051, delimiter=",")
    np.savetxt("s10.csv", s10, delimiter=",")
    np.savetxt("s11.csv", s11, delimiter=",")
    np.savetxt("s20.csv", s20, delimiter=",")
    np.savetxt("s21.csv", s21, delimiter=",")
    np.savetxt("s40.csv", s40, delimiter=",")
    np.savetxt("s41.csv", s41, delimiter=",")
