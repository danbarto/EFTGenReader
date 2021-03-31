import matplotlib.pyplot as plt
import collections
import copy
import numpy as np
import json
import os

### Info about this script ###
#    - This script is very messy, but includes a lot of functions for parsing WCFits
#    - It is designed to take as input text files of WCFits (dumped with the "save" method) and to perform various checks and/or print information about them
#    - Some of the functions are potentially useful in general
#    - Some of the functions are far too specific to be useful in general
#    - So as it stands, this script is probably not really useful, but maybe some variant of it could someday be useful

PROC_NAMES_SHORT = {
    "ttHJet"                  : "ttHJet",
    "ttlnuJet"                : "ttlnuJet",
    "ttllNuNuJetNoHiggs"      : "ttllJet",
    "ttbarJet"                : "ttbarJet",
    "tHq4f"                   : "tHq",
    "tllq4fNoSchanWNoHiggs0p" : "tllq"
}

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

# For run with top19001 lo pt
run5_2heavy2light = {
    "cQq13" : 4.0,
    "cQq83" : 6.0,
    "cQq11" : 23.0,
    "ctq1"  : -19.0,
    "cQq81" : -13.0,
    "ctq8"  : 10.0,
}

# For run with top19001 hi pt
run6_2heavy2light = {
    "cQq13" : -3.0,
    "cQq83" : 4.0,
    "cQq11" : -10.0,
    "ctq1"  : 10.0,
    "cQq81" : 17.0,
    "ctq8"  : -13.0,
}

for19001base_2heavy2light = {
    "cQq13" : 3.0,
    "cQq83" : -4.0,
    "cQq11" : 15.0,
    "ctq1"  : -12.0,
    "cQq81" : 15.0,
    "ctq8"  : -8.0,
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
top19001_ttZJet_pt = {
    "ctW"     : -5.02,
    "ctp"     : 32.91,
    "cpQM"    : -8.06,
    "ctei"    : 6.003,
    "ctli"    : 10.16,
    "cQei"    : 4.804,
    "ctZ"     : -3.86,
    "cQlMi"   : -7.15,
    "cQl3i"   : -8.33,
    "ctG"     : 1.606,
    "ctlTi"   : 2.824,
    "cbW"     : -3.82,
    "cpQ3"    : -13.2,
    "cptb"    : 14.49,
    "cpt"     : -32.6,
    "ctlSi"   : -7.06,
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
    "cQq13" : 0.0,
    "cQq83" : 0.0,
    "cQq11" : 0.0,
    "ctq1"  : 0.0,
    "cQq81" : 0.0,
    "ctq8"  : 0.0,
    "sm" : 1.,
}


WC_LST = [ "ctW" , "ctZ" , "ctp" , "cpQM" , "ctG" , "cbW" , "cpQ3" , "cptb" , "cpt" , "cQl3i" , "cQlMi" , "cQei" , "ctli" , "ctei" , "ctlSi" , "ctlTi" , "cQq13", "cQq83", "cQq11", "ctq1" , "cQq81", "ctq8" ]

### General helper functions ###

# Get the percent difference between two numbers
def get_pdiff(a,b):
    if a==0 and b==0:
        p = 0.0
    elif a+b==0:
        print "WARNING: Trying to find percent difference between two equal and opposite values:",a,b
        p = None
    else:
        #p = 100.*abs((float(a)-float(b))/((float(a)+float(b))/2.0))
        p = 100.*((float(a)-float(b))/((float(a)+float(b))/2.0))
    return p

# Take a wc pt dictionary, retrun a wcpt string
def make_wcpt_str(wcpt_dict):
    s = "rwgt"
    for wc,val in sorted(wcpt_dict.iteritems()):
        if wc == "sm": continue
        s = s + "_" + wc + "_" + str(val)
    return s

# Open a json file, given a name and a path
def open_json(fpath,file_name):
    fname = file_name
    split_str = fname.split(".")
    if split_str[-1] != "json":
        fname = file_name+".json"
    path_to_file = os.path.join(fpath,fname)
    with open (path_to_file,"r") as f:
        d = json.load(f)
    return d

# Combine two dictionaries
def combine_dicts(d1,d2):
    ret_dict = {}
    for k in d1.keys():
        if k in d2:
            print "\nCannot combine these dictionaries, there is an overlap! Exiting...\n"
            raise Exception
    for k,v in d1.iteritems():
        ret_dict[k] = v
    for k,v in d2.iteritems():
        ret_dict[k] = v
    return ret_dict

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
        if wc1 not in wcpt_dict.keys():
            if wc1 == "sm":
                wc1_val = 1.0
            else:
                print "WARNING: No value specified for WC {wc}. Setting it to 0.".format(wc=wc1)
                wc1_val = 0.0
        else:
            wc1_val = wcpt_dict[wc1]
        if wc2 not in wcpt_dict.keys():
            if wc2 == "sm":
                wc2_val = 1.0
            else:
                print "WARNING: No value specified for WC {wc}. Setting it to 0.".format(wc=wc2)
                wc2_val = 0.0
        else:
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

# Scale a fit to its SM value, return the scaled fit
def scale_fit_to_sm(f):
    sm_val = f["sm*sm"]
    return scale_fit(f,1.0/sm_val)

# Scale a fit by a number, return the scaled fit
def scale_fit(f,scale):
    scale_fit = {}
    for term in f.keys():
        scale_fit[term] = f[term]*float(scale)
    return scale_fit

# Takes a multi-dim fit, and returns just the 1d quadratic fit coeffs for one WC
def get_1d_fit(fit_dict,wc):
    s0 = fit_dict["sm*sm"]
    s1 = fit_dict[wc+"*sm"]
    s2 = fit_dict[wc+"*"+wc]
    return [s0,s1,s2]

# Evaluate a 1d quadratic at a given point
def eval_1d_quad(quad_params,x):
    if len(quad_params) != 3:
        print "\nYou are using this function wrong. There should be three params in a quad fit. Exiting...\n"
        raise Exception
    y = float(quad_params[0]) + float(quad_params[1])*float(x) + float(quad_params[2])*float(x)*float(x)
    return y

# Takes as input 1d quadratic fit params, and returns the x value where y crosses some threshold
def find_where_fit_crosses_threshold(quad_params,threshold):
    threshold = float(threshold)
    y = 0
    x = 0
    for x in np.linspace(0,100,10001):
        y_p = eval_1d_quad(quad_params,x)
        y_n = eval_1d_quad(quad_params,-x)
        if y_p >= threshold:
            y = y_p
            break
        if y_n >= threshold:
            y = y_n
            x = -x
            break
    return x,y

# Gives you the name for a summary tree reader sample given proc, tag, run
def reconstruct_sample_name(p,c,r):
    sample_name = "output_{p}_{c}_{r}".format(p=p,c=c,r=r)
    return sample_name


### More complex helper functions, probably only useful in specific cases ###

# Takes as input a path to a dir (that should contain dumped fits), and then loops over all files in the dir and puts into dict
def put_all_fits_from_a_dir_into_dict(dir_path):
    d = {}
    for f in os.listdir(dir_path):
        p,c,r = (f.split(".")[0]).split("_") # Split the name (not including .txt extension) into process, coeff, run
        #tag = p+"_"+r[-1] # e.g. ttH_test_run0 -> ttH_0
        tag = p+"_"+r # e.g. ttH_test_run0 -> ttH_0
        d[tag] = read_fit_file_dump_dict(os.path.join(dir_path,f))
    return d

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
            if p is not None:
                if ((abs(p)>threshold) and ((abs(v1)>5e-4) or (abs(v2)>5e-4)) ):
                    print "\t",get_string_for_comp_2_vals(t,v1,v2)
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

# Prints some info about weights at some point
# TODO: should use get_string_for_comp_2_vals
def print_fit_eval_info(tag,orig_wgt,rwgt_wgt):
    pdiff = get_pdiff(rwgt_wgt,orig_wgt)
    s = "\t{tag}: orig wgt, rwgt wgt: {orig_wgt} {rwgt_wgt} -> {p}".format(tag=tag,orig_wgt=orig_wgt,rwgt_wgt=rwgt_wgt,p=pdiff)
    print s

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
    print ""
    #small = -1.
    small = 0.0001
    #f1 = scale_fit(f1,1.0/f1["sm*sm"])
    #f2 = scale_fit(f2,1.0/f2["sm*sm"])
    for k in f1.keys():
        v1 = f1[k]
        v2 = f2[k]
        if abs(v1)>small or abs(v2)>small:
            p = get_pdiff(v1,v2)
            #if abs(p) > -1:
            if abs(p) > 5:
                print "\t",k,"\t1,2+3:",v1,",",v2,"\t->",p
    print ""

# Takes two fits, evaluates them at a particular point, and then prints the cross sections and the percent diff between them
def eval_two_fits_and_print_pdiff(f1,f2,point):
    x1 = eval_fit(f1,point)
    x2 = eval_fit(f2,point)
    p = get_pdiff(x1,x2)
    print "\txsec1,xsec2+3:",x1,",",x2,"\t->",p

# Takes a dictionary of fit dictionaires, loops through WCs, finds where 1d fits pass some threshold, dumps to json
def find_start_pt_and_dump_to_JSON(fit_dicts,wc_lst,threshold,tag):
    start_pts_dict = {}
    for sample_name , fit_dict in fit_dicts.iteritems():
        start_pts_dict[sample_name] = {}
        for wc in wc_lst:
            quad1d = get_1d_fit(fit_dict,wc)
            start_pts_dict[sample_name][wc] = round(find_where_fit_crosses_threshold(quad1d,threshold)[0],2)
    print "\n", tag, "\n" , start_pts_dict
    with open(tag+".json","w") as outfile:
        json.dump(start_pts_dict, outfile)

# Compare rwgt points to orig points
# For the new start pt scan fits, and the "base" fits (i.e. TOP-19-001 start pts)
def do_start_pt_comp(proc_lst, run_lst, all_start_pts, fit_dict_scan, orig_wgts_scan, sample_tag_scan, fit_dict_base, orig_wgts_base, sample_tag_base):

    # Evaluate the fits at starting points for each sample
    for p in proc_lst:
        print "\n",p,"\n"

        # Get fit and orig weight for pt based on TOP-19-001 start
        fit_base = fit_dict_base[p+"_0"]
        sname_base = reconstruct_sample_name(p,sample_tag_base,"run0")
        orig_wgt_base = orig_wgts_base[sname_base]["xsecAtStartScaleToSM"]

        # Loop over the 7 new start pts
        for run in run_lst:
            print run
            # Check all22WCsStartPtCheckV2dim6TopMay20GST samples
            sname = reconstruct_sample_name(p,sample_tag_scan,run)
            #orig_wgt = orig_wgts_scan[sname] # If using the manually produced json
            orig_wgt = orig_wgts_scan[sname]["xsecAtStartScaleToSM"] # For the automatically produced json
            for tag,fit in sorted(fit_dict_scan.iteritems()):
                if p not in tag: continue
                rwgt_wgt = eval_fit(fit,all_start_pts[run][PROC_NAMES_SHORT[p]])
                print_fit_eval_info(tag,orig_wgt,rwgt_wgt)
            # Check all22WCsBaselineStartPtTOP19001dim6TopMay20GST samples
            rwgt_wgt_19001base = eval_fit(fit_base,all_start_pts[run][PROC_NAMES_SHORT[p]])
            print_fit_eval_info(p+" base",orig_wgt,rwgt_wgt_19001base)

        # Also include  the start pt based on TOP-19-001 start pt
        print "top19001base"
        for tag,fit in fit_dict_scan.iteritems():
            if p not in tag: continue
            rwgt_wgt = eval_fit(fit,pt_top19001base[PROC_NAMES_SHORT[p]])
            print_fit_eval_info(tag,orig_wgt_base,rwgt_wgt)
        rwgt_wgt_19001base = eval_fit(fit_base,pt_top19001base[PROC_NAMES_SHORT[p]])
        print_fit_eval_info(p+" base",orig_wgt_base,rwgt_wgt_19001base)

# Takes a tag and two values, returns string: "tag: v1 , v1 -> pdiff(v1,v2)"
def get_string_for_comp_2_vals(tag,v1,v2,threshold=None):
    p = get_pdiff(v1,v2)
    s = "{t}: {v1} , {v2} -> {p}".format(t=tag,v1=v1,v2=v2,p=p)
    if threshold is not None:
        if abs(p) > threshold:
            s = s + "\t\t<- NOTE perc diff greater than " + str(threshold)
    return s


# Compare two sets of fits a list of comparison points
def comp_fits_at_pts(fits1,fits2,comp_pts):
    sample_lst = []
    for sample_name1 in fits1.keys():
        if sample_name1 in fits2:
            sample_lst.append(sample_name1)
    for s in sorted(sample_lst):
        proc,run = s.split("_")
        print s
        #print fits1[s],"\n"
        fits1_sm = scale_fit_to_sm(fits1[s])
        fits2_sm = scale_fit_to_sm(fits2[s])
        for run,point_dict in comp_pts.iteritems():
            # TODO: Could probably use eval_two_fits_and_print_pdiff
            wcpt = comp_pts[run][PROC_NAMES_SHORT[proc]]
            xsec1 = eval_fit(fits1[s],wcpt)
            xsec2 = eval_fit(fits2[s],wcpt)
            print "\t"+get_string_for_comp_2_vals(run,xsec1,xsec2)

# Takes two info dictionaries (from the info in the jsons dumped by runGridpackValidation) and finds common processes, runs
def find_common_keys_in_fit_info_dicts(info_dict1,info_dict2,tag1,tag2,ref_run=None):
    # If using a ref run, then want to compare all fits in 2 to just the ref fits in 1
    proc_lst = []
    run_lst = []
    for k2 in info_dict2.keys():
        if "timestamp" in k2: continue
        p = k2.split("_")[1]
        r = k2.split("_")[3]
        if ref_run is not None:
            k1 = reconstruct_sample_name(p,tag1,ref_run)
        else:
            k1 = reconstruct_sample_name(p,tag1,r)
        # Append to the proc lst and run lst
        if k1 in info_dict1.keys():
            if p not in proc_lst:
                proc_lst.append(p)
            if r not in run_lst:
                run_lst.append(r)
    return sorted(proc_lst),sorted(run_lst)

# Compare info about one sample 1 to a ref from sample 1 (from the info in the jsons dumped by runGridpackValidation)
def comp_sample_info_json_to_ref(info_dict1,info_dict2,tag1,tag2,ref_run):
    #comp_var = "fitEvalAtSMScaleToLSs"
    comp_var = "fitEvalAtSM"
    proc_lst,run_lst = find_common_keys_in_fit_info_dicts(info_dict1,info_dict2,tag1,tag2,ref_run=ref_run)
    print "\nComparing samples with tags:",tag1,tag2,"\n"
    # Loop through processes and do the comparison for each run in sample 2 (against the ref sample from sample 1)
    for p in proc_lst:
        print p
        for r in sorted(run_lst):
            k1 = reconstruct_sample_name(p,tag1,ref_run)
            k2 = reconstruct_sample_name(p,tag2,r)
            if info_dict1[k1]["nLSs"] != info_dict2[k2]["nLSs"]:
                print "ERROR: Trying to compare two samples that do not have the same number of LSs."
                raise Exception
            fit_at_sm_scale_to_ls_1 = info_dict1[k1][comp_var]
            fit_at_sm_scale_to_ls_2 = info_dict2[k2][comp_var]
            print "\t{r} {info}".format(r=r,info=get_string_for_comp_2_vals("fit at sm",fit_at_sm_scale_to_ls_1,fit_at_sm_scale_to_ls_2,5))


# Compare info about two samples (from the info in the jsons dumped by runGridpackValidation)
def comp_sample_info_jsons(info_dict1,info_dict2,tag1,tag2):
    print "\nComparing samples with tags:",tag1,tag2,"\n"
    proc_lst,run_lst = find_common_keys_in_fit_info_dicts(info_dict1,info_dict2,tag1,tag2)
    # Compare the info in the dicts
    for p in proc_lst:
        print p
        for r in run_lst:
            k1 = reconstruct_sample_name(p,tag1,r)
            k2 = reconstruct_sample_name(p,tag2,r)
            if info_dict1[k1]["nLSs"] != info_dict2[k2]["nLSs"]:
                print "ERROR: Trying to compare two samples that do not have the same number of LSs."
                raise Exception
        for r in run_lst:
            k1 = reconstruct_sample_name(p,tag1,r)
            k2 = reconstruct_sample_name(p,tag2,r)
            orig_wgt_1 = info_dict1[k1]["originalXWGTUP"]
            orig_wgt_2 = info_dict2[k2]["originalXWGTUP"]
            #print "\t{r} {info}".format(r=r,info=get_string_for_comp_2_vals("orig wgt",orig_wgt_1,orig_wgt_2,3))
        for r in run_lst:
            k1 = reconstruct_sample_name(p,tag1,r)
            k2 = reconstruct_sample_name(p,tag2,r)
            fit_at_sm_1 = info_dict1[k1]["fitEvalAtSM"]
            fit_at_sm_2 = info_dict2[k2]["fitEvalAtSM"]
            print "\t{r} {info}".format(r=r,info=get_string_for_comp_2_vals("fit at sm",fit_at_sm_1,fit_at_sm_2,10))
        for r in run_lst:
            k1 = reconstruct_sample_name(p,tag1,r)
            k2 = reconstruct_sample_name(p,tag2,r)
            fit_at_sm_scale_to_ls_1 = info_dict1[k1]["fitEvalAtSMScaleToLSs"]
            fit_at_sm_scale_to_ls_2 = info_dict2[k2]["fitEvalAtSMScaleToLSs"]
            #print "\t{r} {info}".format(r=r,info=get_string_for_comp_2_vals("fit at sm scale to ls",fit_at_sm_scale_to_ls_1,fit_at_sm_scale_to_ls_2))
        for r in run_lst:
            k1 = reconstruct_sample_name(p,tag1,r)
            k2 = reconstruct_sample_name(p,tag2,r)
            xsec_at_start_scale_to_ls_1 = info_dict1[k1]["xsecAtStartScaleToLSs"]
            xsec_at_start_scale_to_ls_2 = info_dict2[k2]["xsecAtStartScaleToLSs"]
            #print "\t{r} {info}".format(r=r,info=get_string_for_comp_2_vals("xsec at start scale to ls",xsec_at_start_scale_to_ls_1,xsec_at_start_scale_to_ls_2))

pt_top19001base = {
    "ttHJet"   : combine_dicts(top19001_ttHJet_pt,for19001base_2heavy2light),
    "ttlnuJet" : combine_dicts(top19001_ttWJet_pt,for19001base_2heavy2light),
    "ttllJet"  : combine_dicts(top19001_ttZJet_pt,for19001base_2heavy2light),
    "ttbarJet" : combine_dicts(top19001_ttZJet_pt,for19001base_2heavy2light),
    "tllq"     : combine_dicts(top19001_ttHJet_pt,for19001base_2heavy2light),
    "tHq"      : combine_dicts(top19001_ttZJet_pt,for19001base_2heavy2light),
}

### Main functions ###

# This was the main function for what this script was originally designed for (i.e. trying to understand differences between fits from dim6=1 and dim6^2=2 samples)
def main_for_dim6SQ_checks():

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

    # This should be cleaned up...
    file_names = {
        ### Remake ttHtoLL sample with mmll=10 to be consistnet with the ttll samples
        # Should be from ttHTOll_testHiggsInterferenceWithllWithMMLL10dim6TopMay20GST_GEN_UL17-GEN summay tree sample
        "ttHtollMLL10" : "fit_coeffs/H_int_studies/fitparams_testttHtollMLL10_noNormScalby500_ttHTOllUL17testHiggsInterferenceWithllWithMMLL10dim6TopMay20GSTMatchOffrun0.txt",
        ### Samples used for studies shown in March 5 2021 EFT meeting ###
        # Should be from ttX-tXq_dim6TopMay20GST_testHiggsInterferenceWithLL_GEN_UL17-GEN summary tree sample
        "ttll"    : "fit_coeffs/H_int_studies/fitparams_test_noNormScalby500_ttllUL17testHiggsInterferenceWithlldim6TopMay20GSTrun0.txt",
        "ttllNoH" : "fit_coeffs/H_int_studies/fitparams_test_noNormScalby500_ttllNoHiggsUL17testHiggsInterferenceWithlldim6TopMay20GSTrun0.txt",
        "ttHtoll" : "fit_coeffs/H_int_studies/fitparams_test_noNormScalby500_ttHTOllUL17testHiggsInterferenceWithlldim6TopMay20GSTMatchOffrun0.txt",
        "tllq"    : "fit_coeffs/H_int_studies/fitparams_test_noNormScalby500_tllq4fNoSchanWUL17testHiggsInterferenceWithlldim6TopMay20GSTrun1.txt",
        "tllqNoH" : "fit_coeffs/H_int_studies/fitparams_test_noNormScalby500_tllq4fNoSchanWNoHiggs0pUL17testHiggsInterferenceWithlldim6TopMay20GSTrun1.txt",
        "tHtollq" : "fit_coeffs/H_int_studies/fitparams_test_noNormScalby500_tHTOllq4fNoSchanWUL17testHiggsInterferenceWithlldim6TopMay20GSTrun0.txt",
    }

    all_fits = {}
    for tag,fname in file_names.iteritems():
        f_dict = read_fit_file_dump_dict(fname)
        all_fits[tag] = f_dict

    pts_lst = [top19001hi_pt,top19001lo_pt,OLD_AN_PT,top19001_ttHJet_pt,top19001_ttWJet_pt,all2_pt,just_ctG2_pt,just_ctGn2_pt,some_problematic_looking_coeffs_on_pt,sm_pt]
    samples_lst = ["ttll","ttllNoH","ttHtoll","tllq","tllqNoH","tHtollq"]
    run_lst = ["","_r0","_r1","_r2"]

    '''
    # Compare the sample with mmll=10 in run card against the one original (which accidnetly had no mmll cut since was using ttH run card)
    print_pdiff_if_one_is_not_small(all_fits["ttHtollMLL10"],all_fits["ttHtoll"])
    '''

    '''
    # Check original samples against new samples at different and same starting points (to see variations are comparable to or larger than original 1 vs 2+3 differences)
    for r in run_lst:
        print r
        for s in samples_lst:
            print s
            print_pdiff_if_one_is_not_small(all_fits[s+r],all_fits[s+"_sep"+r]) # Verify that the samples ran through the plotting code together and seperately produce same restults (they do)
            #print_pdiff_if_one_is_not_small(all_fits[s],all_fits[s+"_sep"+r])
            #print_pdiff_if_one_is_not_small(all_fits[s],all_fits[s+"_tog"+r])
            #print_pdiff_if_one_is_not_small(all_fits[s],all_fits[s+"_togRerun"+r])
    '''

    #'''
    # Check 1 vs 2+3 for all samples:
    for r in run_lst:
        print r
        # Uncomment one ttX1 line and one tXq1 line
        #ttX1 = all_fits["ttll"] ; ttX2p3 = add_fits(all_fits["ttllNoH"],all_fits["ttHtollMLL10"]) # As shown at March 5 2021 group meeting, except with mmll=10 sample
        ttX1 = all_fits["ttll"] ; ttX2p3 = add_fits(all_fits["ttllNoH"],all_fits["ttHtoll"]) # As shown at March 5 2021 group meeting
        tXq1 = all_fits["tllq"] ; tXq2p3 = add_fits(all_fits["tllqNoH"],all_fits["tHtollq"]) # As shown at March 5 2021 group meeting
        #ttX1 = all_fits["ttll_tog"+r] ; ttX2p3 = add_fits(all_fits["ttllNoH_tog"+r],all_fits["ttHtoll_tog"+r])
        #tXq1 = all_fits["tllq_tog"+r] ; tXq2p3 = add_fits(all_fits["tllqNoH_tog"+r],all_fits["tHtollq_tog"+r])
        #ttX1 = all_fits["ttll_togRerun"+r] ; ttX2p3 = add_fits(all_fits["ttllNoH_togRerun"+r],all_fits["ttHtoll_togRerun"+r])
        #tXq1 = all_fits["tllq_togRerun"+r] ; tXq2p3 = add_fits(all_fits["tllqNoH_togRerun"+r],all_fits["tHtollq_togRerun"+r])
        # Compare the fit coeffecients
        print_pdiff_if_one_is_not_small(ttX1,ttX2p3)
        print_pdiff_if_one_is_not_small(tXq1,tXq2p3)
        # Compare the xsec at several points
        for p in pts_lst:
            eval_two_fits_and_print_pdiff(ttX1,ttX2p3,p)
            eval_two_fits_and_print_pdiff(tXq1,tXq2p3,p)
    #'''

# Find the place where the WCs scale each process by a certain amount, to be used as the start point
def main_for_finding_start_pts():
    file_names = {
        "ttHJet"   : "fit_coeffs/all22WCs/fitparams_smNorm_ttHJetUL17all22WCsBaselineStartPtTOP19001dim6TopMay20GSTrun0.txt",
        "ttlnuJet" : "fit_coeffs/all22WCs/fitparams_smNorm_ttlnuJetUL17all22WCsBaselineStartPtTOP19001dim6TopMay20GSTrun0.txt",
        "ttllJet"  : "fit_coeffs/all22WCs/fitparams_smNorm_ttllNuNuJetNoHiggsUL17all22WCsBaselineStartPtTOP19001dim6TopMay20GSTrun0.txt",
        "tllq"     : "fit_coeffs/all22WCs/fitparams_smNorm_tllq4fNoSchanWNoHiggs0pUL17all22WCsBaselineStartPtTOP19001dim6TopMay20GSTrun0.txt",
        "tHq"      : "fit_coeffs/all22WCs/fitparams_smNorm_tHq4fUL17all22WCsBaselineStartPtTOP19001dim6TopMay20GSTrun0.txt",
        "ttbarJet" : "fit_coeffs/all22WCs/fitparams_smNorm_ttbarJetUL17all22WCsBaselineStartPtTOP19001dim6TopMay20GSTrun0.txt",
    }
    all_fits = {}
    for tag,fname in file_names.iteritems():
        f_dict = read_fit_file_dump_dict(fname)
        all_fits[tag] = f_dict

    find_start_pt_and_dump_to_JSON(all_fits,WC_LST,5.0,"startpts_scale_by_5p0")
    find_start_pt_and_dump_to_JSON(all_fits,WC_LST,2.0,"startpts_scale_by_2p0")
    find_start_pt_and_dump_to_JSON(all_fits,WC_LST,1.5,"startpts_scale_by_1p5")
    find_start_pt_and_dump_to_JSON(all_fits,WC_LST,1.3,"startpts_scale_by_1p3")
    find_start_pt_and_dump_to_JSON(all_fits,WC_LST,1.1,"startpts_scale_by_1p1")





# Perform validation checks for the start points and reweighting
def main_for_rwgt_validaiton():

    run_lst = ["run0","run1","run2","run3","run4","run5","run6"]
    p_lst_no_ttll = ["ttHJet","ttlnuJet","tHq4f","tllq4fNoSchanWNoHiggs0p","ttbarJet"]
    p_lst = ["ttHJet","ttlnuJet","ttllNuNuJetNoHiggs","tHq4f","tllq4fNoSchanWNoHiggs0p","ttbarJet"]

    #orig_wgts                      = open_json("fit_coeffs/start_pt_checks/","xsec_at_start_pts_by-hand-cp")
    #orig_wgts_19001base            = open_json("fit_coeffs/start_pt_checks/sample_info/old/","fit_info_all22WCsBaselineStartPtTOP19001")
    #orig_wgts_all22WCsStartPtCheck = open_json("fit_coeffs/start_pt_checks/sample_info/old/","fit_info_all22WCsStartPtCheck")
    #orig_wgts_19001samples         = open_json("fit_coeffs/start_pt_checks/sample_info/","fit_info_top19001_samples")
    info_all22WCsBaslineStart   = open_json("fit_coeffs/start_pt_checks/sample_info/","fit_info_all22WCsBaselineStartPtTOP19001.json")
    info_all22WCsStartPtCheck   = open_json("fit_coeffs/start_pt_checks/sample_info/","fit_info_all22WCsStartPtCheck.json")
    info_all22WCsStartPtCheckV2 = open_json("fit_coeffs/start_pt_checks/sample_info/","fit_info_all22WCsStartPtCheckV2.json")
    info_dim6eq0                = open_json("fit_coeffs/start_pt_checks/sample_info/","fit_info_dim6eq0.json")

    # Put all of the start pts into dictionaries
    all_start_pts = {}
    start_pts_file_lst = ["startpts_scale_by_1p1","startpts_scale_by_1p3","startpts_scale_by_1p5","startpts_scale_by_2p0","startpts_scale_by_5p0"]
    for idx,f_name in enumerate(start_pts_file_lst):
        all_start_pts["run"+str(idx)] = open_json("fit_coeffs/start_pt_checks/start_pt_jsons",f_name)
    # Put in the other two WCs that do not come from jsons (note these are used for all processes)
    all_start_pts["run5"] = {}
    all_start_pts["run6"] = {}
    for p in p_lst:
        all_start_pts["run5"][PROC_NAMES_SHORT[p]] = combine_dicts(top19001lo_pt,run5_2heavy2light)
        all_start_pts["run6"][PROC_NAMES_SHORT[p]] = combine_dicts(top19001hi_pt,run6_2heavy2light)

    # Put the fits into a dictionary
    startPtScan_fits   = put_all_fits_from_a_dir_into_dict("fit_coeffs/start_pt_checks/normToSM/all22WCsStartPtCheck")
    startPtScanV2_fits = put_all_fits_from_a_dir_into_dict("fit_coeffs/start_pt_checks/normToSM/all22WCsStartPtCheckV2")
    startPtBase_fits   = put_all_fits_from_a_dir_into_dict("fit_coeffs/start_pt_checks/normToSM/all22WCsBaselineStartPtTOP19001")
    top19001_fits      = put_all_fits_from_a_dir_into_dict("fit_coeffs/start_pt_checks/normToSM/Top19001Samples")
    fits_all22WCsBaslineStart   = put_all_fits_from_a_dir_into_dict("fit_coeffs/start_pt_checks/normToLS/all22WCsBaselineStartPtTOP19001")
    fits_all22WCsStartPtCheck   = put_all_fits_from_a_dir_into_dict("fit_coeffs/start_pt_checks/normToLS/all22WCsStartPtCheck")
    fits_all22WCsStartPtCheckV2 = put_all_fits_from_a_dir_into_dict("fit_coeffs/start_pt_checks/normToLS/all22WCsStartPtCheckV2")
    fits_dim6eq0                = put_all_fits_from_a_dir_into_dict("fit_coeffs/start_pt_checks/normToLS/dim6Eq0")

    # Compare the v0 and v1 sample info to each other
    comp_sample_info_jsons(info_all22WCsStartPtCheck,info_all22WCsStartPtCheckV2,"all22WCsStartPtCheckdim6TopMay20GST","all22WCsStartPtCheckV2dim6TopMay20GST")

    # Compare v0 and v1 sample info to the dim6eq0 samples
    #comp_sample_info_json_to_ref(info_dim6eq0,info_all22WCsStartPtCheck,"all22WCsDim6Eq0dim6TopMay20GST","all22WCsStartPtCheckdim6TopMay20GST","run0")
    #comp_sample_info_json_to_ref(info_dim6eq0,info_all22WCsStartPtCheckV2,"all22WCsDim6Eq0dim6TopMay20GST","all22WCsStartPtCheckV2dim6TopMay20GST","run0")

    # Compare the v0 and v1 fits to each other
    #comp_fits_at_pts(fits_all22WCsStartPtCheck,fits_all22WCsStartPtCheckV2,all_start_pts)

    # Compare the v0 and v1 coeffecients from LS norm to SM norm
    for sample in startPtScan_fits.keys():
        #print "\n",sample,"\n"
        sm_fit = startPtScan_fits[sample]
        ls_fit = scale_fit_to_sm(fits_all22WCsStartPtCheck[sample])
        #get_fit_dicts_pdiff(sm_fit,ls_fit,10)
        for run,point_dict in all_start_pts.iteritems():
            wcpt = point_dict[PROC_NAMES_SHORT[p]]
            #eval_two_fits_and_print_pdiff(sm_fit,ls_fit,wcpt)


    '''
    ########################################
    ### Compare the v0 to the v1 samples ###

    # Compare the v0 to the v1 samples original weights
    for sample,info in orig_wgts_all22WCsStartPtCheck.iteritems():
        if "timestamp" in sample: continue
        if "ttll" in sample: continue
        _,p,c,r = sample.split("_")
        v2_s_name = reconstruct_sample_name(p,"all22WCsStartPtCheckV2dim6TopMay20GST",r)
        #v0_wgt_orig = info["xsecAtStartScaleToSM"]
        v0_wgt_orig = info["originalXWGTUP"]
        v1_wgt_orig = orig_wgts[v2_s_name]
        p = get_pdiff(v0_wgt_orig,v1_wgt_orig)
        if abs(p) > 2.0:
            print "{s}\n\tv0,v1: {w0},{w1} -> {p}\n".format(s=sample,w0=v0_wgt_orig,w1=v1_wgt_orig,p=p)

    raise Exception
    '''
    '''
    # Compare the samples' reweighted values to each other
    for sample in startPtScanV2_fits.keys():
        print sample
        fit_v0 = startPtScan_fits[sample]
        fit_v1 = startPtScanV2_fits[sample]
        p = sample.split("_")[0]
        r = "run"+sample.split("_")[1]
        for run,point_dict in all_start_pts.iteritems():
            wcpt = point_dict[PROC_NAMES_SHORT[p]]
            wgt0 = eval_fit(fit_v0,wcpt)
            wgt1 = eval_fit(fit_v1,wcpt)
            print "\t",run,wgt0,wgt1,"->",get_pdiff(wgt0,wgt1)

    raise Exception
    ########################################
    '''

    ## Check top 19001 samples
    #for sample,info in orig_wgts_19001samples.iteritems():
    #    if "timestamp" in sample: continue
    #    _,p,c,r = sample.split("_")
    #    orig_wgt = info["xsecAtStartScaleToSM"]
    #    start_pt = pt_top19001base[PROC_NAMES_SHORT[p]]
    #    for sample_name,fit in top19001_fits.iteritems():
    #        rwgt_wgt = eval_fit(fit,start_pt)
    #        print sample_name,orig_wgt,rwgt_wgt

    # Compare the rwgt points to the orig pts
    #do_start_pt_comp(p_lst_no_ttll, run_lst, all_start_pts, startPtScanV2_fits, orig_wgts, "all22WCsStartPtCheckV2dim6TopMay20GST", startPtBase_fits, orig_wgts_19001base, "all22WCsBaselineStartPtTOP19001dim6TopMay20GST")
    #do_start_pt_comp(p_lst, run_lst, all_start_pts, startPtScan_fits, orig_wgts_all22WCsStartPtCheck, "all22WCsStartPtCheckdim6TopMay20GST", startPtBase_fits, orig_wgts_19001base, "all22WCsBaselineStartPtTOP19001dim6TopMay20GST")



# Run one of the main functions
#main_for_dim6SQ_checks()
#main_for_hll_int_checks()
#main_for_finding_start_pts()
main_for_rwgt_validaiton()
