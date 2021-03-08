import matplotlib.pyplot as plt
import collections
import copy

### Info about this script ###
#    - This script is very messy, but includes a lot of functions for parsing WCFits
#    - It is designed to take as input text files of WCFits (dumped with the "save" method) and to perform various checks and/or print information about them
#    - Some of the functions are potentially useful in general
#    - Some of the functions are far too specific to be useful in general
#    - So as it stands, this script is probably not really useful, but maybe some variant of it could someday be useful


### Some WCPoints  ###

top19001hi_pt = {
    "ctW"   : 2.87,
    "ctZ"   : 3.15,
    "ctp"   :  44.26,
    "cpQM"  : 21.65,
    "ctG"   : 1.18,
    "cbW"   : 4.95,
    "cpQ3"  : 3.48,
    "cptb"  :  12.63,
    "cpt"   :  12.31,
    "cQl3i" : 8.97,
    "cQlMi" : 4.99,
    "cQei"  : 4.59,
    "ctli"  : 4.82,
    "ctei"  : 4.86,
    "ctlSi" : 6.52,
    "ctlTi" : 0.84,
    "sm" : 1.,
}
top19001lo_pt = {
    "ctW"   : -3.08,
    "ctZ"   : -3.32,
    "ctp"   : -16.98,
    "cpQM"  : -7.59,
    "ctG"   : -1.38,
    "cbW"   : -4.95,
    "cpQ3"  : -7.37,
    "cptb"  : -12.72,
    "cpt"   : -18.62,
    "cQl3i" : -9.67,
    "cQlMi" : -4.02,
    "cQei"  : -4.38,
    "ctli"  : -4.29,
    "ctei"  : -4.24,
    "ctlSi" : -6.52,
    "ctlTi" : -0.84,
    "sm" : 1.,
}

OLD_AN_PT = {
    "ctW"   : -0.58,
    "ctZ"   : -0.63,
    "ctp"   : 25.50,
    "cpQM"  : -1.07,
    "ctG"   : -0.85,
    "cbW"   : 3.17 ,
    "cpQ3"  : -1.81,
    "cptb"  : 0.13 ,
    "cpt"   : -3.25,
    "cQl3i" : -4.20,
    "cQlMi" : 0.74 ,
    "cQei"  : -0.27,
    "ctli"  : 0.33 ,
    "ctei"  : 0.33 ,
    "ctlSi" : -0.07,
    "ctlTi" : -0.01,
    "sm" : 1.,
}
top19001_ttHJet_pt = {
    "ctW"    :  -8.30  ,
    "ctp"    :  64.33  ,
    "cpQM"   :  45.88  ,
    "ctei"   :  24.32  ,
    "ctli"   :  24.43  ,
    "cQei"   :  23.75  ,
    "ctZ"    :  -6.09  ,
    "cQlMi"  :  23.95  ,
    "cQl3i"  :  21.54  ,
    "ctG"    :  -3.60  ,
    "ctlTi"  :  21.80  ,
    "cbW"    :  49.59  ,
    "cpQ3"   :  -51.1  ,
    "cptb"   :  136.1  ,
    "cpt"    :  -43.5  ,
    "ctlSi"  :  -20.0  ,
    "sm" : 1.,
}

top19001_ttWJet_pt = {
    "ctW"    : -3.82,
    "ctp"    : 51.50,
    "cpQM"   : 23.00,
    "ctei"   : 8.938,
    "ctli"   : -7.00,
    "cQei"   : 8.968,
    "ctZ"    : 5.727,
    "cQlMi"  : 6.952,
    "cQl3i"  : 9.243,
    "ctG"    : 2.430,
    "ctlTi"  : 2.116,
    "cbW"    : -7.37,
    "cpQ3"   : -14.4,
    "cptb"   : -21.8,
    "cpt"    : -20.3,
    "ctlSi"  : -9.99,
    "sm" : 1.,
}
all2_pt = {
    "ctW"   : 2.0,
    "ctZ"   : 2.0,
    "ctp"   : 2.0,
    "cpQM"  : 2.0,
    "ctG"   : 2.0,
    "cbW"   : 2.0,
    "cpQ3"  : 2.0,
    "cptb"  : 2.0,
    "cpt"   : 2.0,
    "cQl3i" : 2.0,
    "cQlMi" : 2.0,
    "cQei"  : 2.0,
    "ctli"  : 2.0,
    "ctei"  : 2.0,
    "ctlSi" : 2.0,
    "ctlTi" : 2.0,
    "sm" : 1.,
}
just_ctG2_pt = {
    "ctW"   : 0.0,
    "ctZ"   : 0.0,
    "ctp"   : 0.0,
    "cpQM"  : 0.0,
    "ctG"   : 2.0,
    "cbW"   : 0.0,
    "cpQ3"  : 0.0,
    "cptb"  : 0.0,
    "cpt"   : 0.0,
    "cQl3i" : 0.0,
    "cQlMi" : 0.0,
    "cQei"  : 0.0,
    "ctli"  : 0.0,
    "ctei"  : 0.0,
    "ctlSi" : 0.0,
    "ctlTi" : 0.0,
    "sm" : 1.,
}
just_ctGn2_pt = {
    "ctW"   : 0.0,
    "ctZ"   : 0.0,
    "ctp"   : 0.0,
    "cpQM"  : 0.0,
    "ctG"   : -2.0,
    "cbW"   : 0.0,
    "cpQ3"  : 0.0,
    "cptb"  : 0.0,
    "cpt"   : 0.0,
    "cQl3i" : 0.0,
    "cQlMi" : 0.0,
    "cQei"  : 0.0,
    "ctli"  : 0.0,
    "ctei"  : 0.0,
    "ctlSi" : 0.0,
    "ctlTi" : 0.0,
    "sm" : 1.,
}

some_problematic_looking_coeffs_on_pt = {
    "ctW"   : -1.5,
    "ctZ"   : 4.0,
    "ctp"   : 0.0,
    "cpQM"  : 0.0,
    "ctG"   : 1.0,
    "cbW"   : 3.0,
    "cpQ3"  : 5.0,
    "cptb"  : 0.0,
    "cpt"   : 15.0,
    "cQl3i" : 10.0,
    "cQlMi" : 0.0,
    "cQei"  : 5.0,
    "ctli"  : 0.0,
    "ctei"  : 0.0,
    "ctlSi" : -7.0,
    "ctlTi" : 0.0,
    "sm" : 1.,
}

sm_pt = {
    "ctW"   : 0.0,
    "ctZ"   : 0.0,
    "ctp"   : 0.0,
    "cpQM"  : 0.0,
    "ctG"   : 0.0,
    "cbW"   : 0.0,
    "cpQ3"  : 0.0,
    "cptb"  : 0.0,
    "cpt"   : 0.0,
    "cQl3i" : 0.0,
    "cQlMi" : 0.0,
    "cQei"  : 0.0,
    "ctli"  : 0.0,
    "ctei"  : 0.0,
    "ctlSi" : 0.0,
    "ctlTi" : 0.0,
    "sm" : 1.,
}

### General helper functions ###

# Get the percent difference between two numbers
def get_pdiff(a,b):
    if a+b==0:
        p = None
    elif a==0 and b==0:
        p = 0.0
    else:
        #p = 100.*abs((float(a)-float(b))/((float(a)+float(b))/2.0))
        p = 100.*((float(a)-float(b))/((float(a)+float(b))/2.0))
    return p

# Open a text file that we've dumped a WCFit into and put the info into a dictionary
def read_fit_file_dump_dict(fpath):
    with open (fpath) as f:
        f_lst = [] # The terms should be the first entry in this list, the values of those terms are the second
        fit_dict = {}
        for l in f:
            f_lst.append(l)
        coeff_terms_lst = f_lst[0].split()
        coeff_vals_lst = f_lst[1].split()
        for i in range(len(coeff_terms_lst)):
            fit_dict[coeff_terms_lst[i]] = float(coeff_vals_lst[i+1]) # The values list starts with the tag name
        return fit_dict

# Evalueate a fit dictionary at some point in the 16d phase space
def eval_fit(fit_dict,wcpt_dict):
    xsec = 0
    for wc_str, coeff_val in fit_dict.iteritems():
        wc1,wc2 = wc_str.split("*")
        wc1_val = wcpt_dict[wc1]
        wc2_val = wcpt_dict[wc2]
        xsec = xsec + wc1_val*wc2_val*coeff_val
    return xsec

# Adds two fits together, return the sum
def add_fits(f1,f2):
    if collections.Counter(f1.keys()) != collections.Counter(f2.keys()): # Order of keys does not matter
        print "\n Canot compare these fit dictionaries, they do not have the same keys!"
        raise Exception
    sum_fit = {}
    for term in f1.keys():
        #print f1[term] + f2[term]
        sum_fit[term] = f1[term] + f2[term]
    return sum_fit

# Scale a fit by a number, return the scaled fit
def scale_fit(f,scale):
    scale_fit = {}
    for term in f.keys():
        scale_fit[term] = f[term]*float(scale)
    return scale_fit


### More complex helper functions ###

# Takes two fit dictionaries and returns a dictioinary with the percent difference between the values
def get_fit_dicts_pdiff(d1,d2,threshold=None):
    terms_lst1 = d1.keys()
    terms_lst2 = d2.keys()
    if collections.Counter(terms_lst1) != collections.Counter(terms_lst2):
        print "\n Canot compare these dictionaries, they do not have the same keys!"
        raise Exception
    p_dict = {}
    for t in terms_lst1:
        v1 = d1[t]
        v2 = d2[t]
        p = get_pdiff(v1,v2)
        p_dict[t]=p
        if threshold is not None:
            if p>threshold:
                print t,v1,v2,p
        else:
            print t,v1,v2,p
    return p_dict

# Takes two fit dictionaries and fills hists for percent differences between them
def make_hist_of_diffs(d1,d2,savename):
    pdiff_dict = get_fit_dicts_pdiff(d1,d2)
    #print p_dict
    pdiff_lin_lst = []
    pdiff_quad_lst = []
    for tag,p in pdiff_dict.iteritems():
        if "sm" in tag:
            pdiff_lin_lst.append(p)
        else:
            pdiff_quad_lst.append(p)
    print "max, avg quad diff:",max(pdiff_quad_lst),sum(pdiff_quad_lst)/len(pdiff_quad_lst)
    print "max,avg  lin diff:",max(pdiff_lin_lst),sum(pdiff_lin_lst)/len(pdiff_lin_lst)
    # Get rid of 0s from quad:
    pdiff_quad_lst_no0s = []
    for i in pdiff_quad_lst:
        if i>0:
            pdiff_quad_lst_no0s.append(i)
    plt.hist(pdiff_lin_lst,bins=40)
    plt.title(savename)
    plt.xlabel('Percent diff between DIM6=1 and DIM6^2 linear terms')
    plt.savefig(savename+"_lin.png")
    plt.clf()
    #plt.hist(pdiff_quad_lst,bins=40)
    plt.hist(pdiff_quad_lst_no0s,bins=40)
    plt.title(savename)
    plt.xlabel('Percent diff between DIM6=1 and DIM6^2 linear terms')
    plt.xlabel('Percent diff between DIM6=1 and DIM6^2 quad terms')
    plt.savefig(savename+"_quad.png")
    plt.clf()


### Functions that print info about things ###

# Print the values in a fit dictionary
def print_fit_dict(d,threshold=None):
    for term,val in d.iteritems():
        if threshold is not None:
            if abs(val)>threshold:
                print term, val
        else:
            print term,val

# Print the coeffecients for a particular WC in the dictionary
def print_by_wc(wc,d):
    for term,val in d.iteritems():
        quad_str = wc+"*"+wc
        if (quad_str == term):
            print "quad:",term,val
        elif ((wc in term) and ("sm" in term)):
            print "lin:",term,val

# Print a bunch of info about a percent difference dictionary
def print_pdiff_info(p_dict,threshold=None):
    lin_dict = {}
    quad_dict = {}
    # Loop over terms and pdiffs and print values
    print "Pdiff values over",threshold,":"
    for t,p in p_dict.iteritems():
        if threshold is not None:
            if p>threshold:
                print "\t",t,p
        else:
            print "\t",t,p
        if "sm" in t:
            lin_dict[t] = p
        else:
            quad_dict[t] = p
    # Print info about linear terms:
    print "Linear info:"
    lmax = max(lin_dict.values())
    lavg = sum(lin_dict.values())/len(lin_dict.values())
    print "\tmax and ave:",lmax,lavg
    for t,p in lin_dict.iteritems():
        if p !=0:
            print "\tNon zero lin:",t,p
    # Print info about quad terms:
    print "Quad info:"
    qmax = max(quad_dict.values())
    qavg = sum(quad_dict.values())/len(lin_dict.values())
    print "\tmax and ave:",qmax,qavg
    for t,p in quad_dict.iteritems():
        if p !=0:
            print "\tNon zero quad:",t,p


### Very specific functions that are probably not useful for anything besides the very specifc case they were originally used for ###

# Print some info about two fits, and make histograms
def comp_by_tags(tag1,tag2,all_fits_d,threshold):
    print "\nComparing:",tag1,tag2,",threshold=",threshold
    pdiff_dict = get_fit_dicts_pdiff(all_fits_d[tag1],all_fits_d[tag2],threshold)
    print_pdiff_info(pdiff_dict,30)
    print ""
    make_hist_of_diffs(all_fits_d[tag1],all_fits_d[tag2],tag1+tag2)

# Do some checks about what causes the differences between two fits
def explore_effects(fit_dict1,fit_dict2,wcpt_dict):
    wcpt_dict = copy.deepcopy(wcpt_dict)
    fit_dict1 = copy.deepcopy(fit_dict1)
    fit_dict2 = copy.deepcopy(fit_dict2)
    print "Before:"
    f1 = eval_fit(fit_dict1,wcpt_dict)
    f2 = eval_fit(fit_dict2,wcpt_dict)
    print f1,f2,get_pdiff(f1,f2)
    wcpt_dict["ctG"]=0
    wcpt_dict["ctp"]=0
    #top19001hi_pt["ctp"]=0
    #top19001hi_pt["ctZ"]=0
    #top19001hi_pt["ctW"]=0
    #print all_fits["eq1_run0"]["ctp*ctG"]
    #print all_fits["sq_run0"]["ctp*ctG"]
    #fit_dict1["ctp*ctG"] = fit_dict2["ctp*ctG"]
    #tmpf = fit_dict1
    #tmpf["ctp*ctG"] = fit_dict2["ctp*ctG"]
    print "After:"
    f1 = eval_fit(fit_dict1,wcpt_dict)
    f2 = eval_fit(fit_dict2,wcpt_dict)
    print f1,f2,get_pdiff(f1,f2)

# Print the p diff between two fit dictionaries if at least one of them is not tiny
def print_pdiff_if_one_is_not_small(f1,f2):
    small = 0.0001
    for k in f1.keys():
        v1 = f1[k]
        v2 = f2[k]
        #if v1>small or v2>small:
        if abs(v1)>small or abs(v2)>small:
            p = get_pdiff(v1,v2)
            if abs(p) > 5:
                print "\t",k,"\t1,2+3:",v1,",",v2,"\t->",p

# Takes two fits, evaluates them at a particular point, and then prints the cross sections and the percent diff between them
def eval_two_fits_and_print_pdiff(f1,f2,point):
    x1 = eval_fit(f1,point)
    x2 = eval_fit(f2,point)
    p = get_pdiff(x1,x2)
    print "\txsec1,xsec2+3:",x1,",",x2,"\t->",p


### Main functions ###

# This was the main function for what this script was originally designed for (trying to understand differences between fits from dim6=1 and dim6^2=2 samples)
def main():

    file_names = {
        "test" : "test.txt",
        "eq1_run0" : "fitparams_tthUL17testDIM6EQ1dim6TopMay20GSTStartPtChecksrun0.txt",
        "eq1_run1" : "fitparams_tthUL17testDIM6EQ1dim6TopMay20GSTStartPtChecksrun1.txt",
        "eq1_run2" : "fitparams_tthUL17testDIM6EQ1dim6TopMay20GSTStartPtChecksrun2.txt",
        "eq1_run3" : "fitparams_tthUL17testDIM6EQ1dim6TopMay20GSTStartPtChecksrun3.txt",
        "sq_run0" : "fitparams_tthUL17testDIM6SQdim6TopMay20GSTStartPtChecksrun0.txt",
        "sq_run1" : "fitparams_tthUL17testDIM6SQdim6TopMay20GSTStartPtChecksrun1.txt",
        "sq_run3" : "fitparams_tthUL17testDIM6SQdim6TopMay20GSTStartPtChecksrun3.txt",
        "sq_run2" : "fitparams_tthUL17testV2DIM6SQdim6TopMay20GSTStartPtChecksrun2.txt",
    }

    fit_names_ttHJetctGctp = {
        "eq1_run0" : "fitparams_test_ttHJetUL17testV2JustctGctpDIM6EQ1dim6TopMay20GSTrun0.txt",
        "eq1_run1" : "fitparams_test_ttHJetUL17testV2JustctGctpDIM6EQ1dim6TopMay20GSTrun1.txt",
        "sq_run0" : "fitparams_test_ttHJetUL17testV2JustctGctpDIM6SQdim6TopMay20GSTrun0.txt",
        "sq_run1" : "fitparams_test_ttHJetUL17testV2JustctGctpDIM6SQdim6TopMay20GSTrun1.txt",
    }

    all_fits = {}
    for tag,fname in file_names.iteritems():
        f_dict = read_fit_file_dump_dict(fname)
        all_fits[tag] = f_dict

    all_fits_ttHJetctGctp = {}
    for tag,fname in fit_names_ttHJetctGctp.iteritems():
        f_dict_tmp = read_fit_file_dump_dict(fname)
        all_fits_ttHJetctGctp[tag] = f_dict_tmp

    print_fit_dict(all_fits["eq1_run0"],.009)
    print ""
    print "ctZ"
    print_by_wc("ctZ",all_fits["eq1_run0"])
    print "ctW"
    print_by_wc("ctW",all_fits["eq1_run0"])
    print "ctG"
    print_by_wc("ctG",all_fits["eq1_run0"])
    print "ctp"
    print_by_wc("ctp",all_fits["eq1_run0"])

    '''
    # Compare fits and make plots
    #print all_fits
    #comp_by_tags("eq1_run0","sq_run0",all_fits,25)
    #comp_by_tags("eq1_run1","sq_run1",all_fits,25)
    #comp_by_tags("eq1_run2","sq_run2",all_fits,25)
    #comp_by_tags("eq1_run3","sq_run3",all_fits,25)

    explore_effects(all_fits["eq1_run0"],all_fits["sq_run0"],top19001hi_pt)
    print ""
    explore_effects(all_fits["eq1_run0"],all_fits["sq_run0"],OLD_AN_PT)
    print ""
    explore_effects(all_fits["eq1_run0"],all_fits["sq_run0"],top19001_ttHJet_pt)
    print ""
    explore_effects(all_fits["eq1_run0"],all_fits["sq_run0"],top19001_ttWJet_pt)
    '''
    '''
    # Check the sum and difference between ctG ctp fits so that we can plug info wolphram to visualize
    print_fit_dict(all_fits_ttHJetctGctp["eq1_run0"],0)
    f_sum = add_fits(all_fits_ttHJetctGctp["eq1_run0"],all_fits_ttHJetctGctp["sq_run0"])
    f_dif = add_fits(all_fits_ttHJetctGctp["eq1_run0"],scale_fit(all_fits_ttHJetctGctp["sq_run0"],-1))
    print "\nsum:"
    print_fit_dict(f_sum)
    print "\ndif:"
    print_fit_dict(f_dif)
    '''

# This was the main function for this scripts second use case: Trying to understand if we are missing interference by seperating the diagrams with H
def main_for_hll_int_checks():
    file_names = {
        #"ttll"    : "fitparams_test_noNormNoScalby500_ttllUL17testHiggsInterferenceWithlldim6TopMay20GSTrun0.txt",
        #"ttllNoH" : "fitparams_test_noNormNoScalby500_ttllNoHiggsUL17testHiggsInterferenceWithlldim6TopMay20GSTrun0.txt",
        #"ttHtoll" : "fitparams_test_noNormNoScalby500_ttHTOllUL17testHiggsInterferenceWithlldim6TopMay20GSTMatchOffrun0.txt",
        #"tllq"    : "fitparams_test_noNormNoScalby500_tllq4fNoSchanWUL17testHiggsInterferenceWithlldim6TopMay20GSTrun1.txt",
        #"tllqNoH" : "fitparams_test_noNormNoScalby500_tllq4fNoSchanWNoHiggs0pUL17testHiggsInterferenceWithlldim6TopMay20GSTrun1.txt",
        #"tHtollq" : "fitparams_test_noNormNoScalby500_tHTOllq4fNoSchanWUL17testHiggsInterferenceWithlldim6TopMay20GSTrun0.txt",
        "ttll"    : "fit_coeffs/fitparams_test_noNormScalby500_ttllUL17testHiggsInterferenceWithlldim6TopMay20GSTrun0.txt",
        "ttllNoH" : "fit_coeffs/fitparams_test_noNormScalby500_ttllNoHiggsUL17testHiggsInterferenceWithlldim6TopMay20GSTrun0.txt",
        "ttHtoll" : "fit_coeffs/fitparams_test_noNormScalby500_ttHTOllUL17testHiggsInterferenceWithlldim6TopMay20GSTMatchOffrun0.txt",
        "tllq"    : "fit_coeffs/fitparams_test_noNormScalby500_tllq4fNoSchanWUL17testHiggsInterferenceWithlldim6TopMay20GSTrun1.txt",
        "tllqNoH" : "fit_coeffs/fitparams_test_noNormScalby500_tllq4fNoSchanWNoHiggs0pUL17testHiggsInterferenceWithlldim6TopMay20GSTrun1.txt",
        "tHtollq" : "fit_coeffs/fitparams_test_noNormScalby500_tHTOllq4fNoSchanWUL17testHiggsInterferenceWithlldim6TopMay20GSTrun0.txt",
    }

    all_fits = {}
    for tag,fname in file_names.iteritems():
        f_dict = read_fit_file_dump_dict(fname)
        all_fits[tag] = f_dict

    pts_lst = [top19001hi_pt,top19001lo_pt,OLD_AN_PT,top19001_ttHJet_pt,top19001_ttWJet_pt,all2_pt,just_ctG2_pt,just_ctGn2_pt,some_problematic_looking_coeffs_on_pt,sm_pt]

    ttll1 = all_fits["ttll"]
    ttll2plus3 = add_fits(all_fits["ttllNoH"],all_fits["ttHtoll"])
    print ""
    print_pdiff_if_one_is_not_small(ttll1,ttll2plus3)
    print ""
    for p in pts_lst:
        eval_two_fits_and_print_pdiff(ttll1,ttll2plus3,p)
    print ""

    tllq1 = all_fits["tllq"]
    tllq2plus3 = add_fits(all_fits["tllqNoH"],all_fits["tHtollq"])
    print_pdiff_if_one_is_not_small(tllq1,tllq2plus3)
    print ""
    for p in pts_lst:
        eval_two_fits_and_print_pdiff(tllq1,tllq2plus3,p)


# Run one of the main functions
#main()
main_for_hll_int_checks()
