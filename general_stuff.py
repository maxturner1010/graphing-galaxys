import numpy as np, h5py
import matplotlib.pyplot as plt
from scipy.optimize import curve_fit

file_directory = "flags_galaxy_counts_v0306.hdf5"
best_results = {}
data_pointer = ["F090W","F115W","F150W","F200W","F277W","F356W","F444W"]

def get_attributes(directory):
    # print the attributes of the file 
    data = file[directory]
    attrs = data.attrs.keys()
    print(attrs)
    # loop through the attributes and print them
    for key in attrs:
        print(f"{key}: {data.attrs[str(key)]}")
def general_attributes(name):
    # print name of file and description if it is needed
    print(name)
    stuff = file[name]
    if 'Description' in stuff.attrs:
        print(f"dectrption:{stuff.attrs['Description']}")
def genpow_func(x,*args):
    total = 0
    for i in range(0,len(args),2):
        total += args[i] * (x ** args[i+1])
    return total
def genpol_func(x,*args):
    # ploynomial function based on the amount or args
    counter = 0
    total = 0
    for var in args:
        total += var * (x ** (counter))
        counter += 1
    return total
def genlog_func(x,*args):
    counter = 0
    total = 0
    logx = np.log(x)
    for var in args:
        total += var * (logx ** (counter))
        counter += 1
    return total
def genpol_print(*args):
    # print general polynomial function based on amount of args
    string = f"y = "
    counter = 0
    for var in args:
        if int(var) > 0:
            string = string + "+ "
        string += f"{var} * x^{counter} "
        counter += 1
    return string
def genlog_print(*args):
    counter = 0
    func = "y = "
    for var in args:
        if var > 0:
            func += "+"
        func += f"{var} x logx^{counter} "
        counter += 1
    return func
def get_chi2(true_data,modle_data,errors):
    #get kai2 for each data point in np array and the find sum
    kai2 = ((true_data - modle_data)/errors)**2
    return np.sum(kai2)
def find_func(bin_centre,s,errors,n,best_results,pointer):
    # searches through polynomials and returns the one with chi2 close to 1
    best_chi2 = -1
    best_aic = -1
    best_bic = -1
    for i in range(2,13):
        try:
            test = np.ones(i)
            var1 = curve_fit(genpol_func,bin_centre[s],n[s],sigma = errors[s],p0 = test)[0]
            modle_data = genpol_func(bin_centre[s],*var1)
            reduce_chi2,aic,bic = check_module(n[s],modle_data,errors[s],i)
            differance = abs(reduce_chi2 - 1)
            if best_chi2 == -1 or best_chi2 > differance:
                best_chi2 = differance
                best_chi2_data = (var1,genpol_func,reduce_chi2,aic,bic)
            if best_aic == -1 or best_aic > aic:
                best_aic = aic
                best_aic_data = (var1,genpol_func,reduce_chi2,aic,bic)
            if best_bic == -1 or best_bic > bic:
                best_bic = bic
                best_bic_data = (var1,genpol_func,reduce_chi2,aic,bic)
        except Exception as e:
            print(e)
            print(i)
    for i in range(1,13):
        try:
            test = np.ones(i)
            var1 = curve_fit(genlog_func,bin_centre[s],n[s],sigma = errors[s],p0 = test)[0]
            modle_data = genlog_func(bin_centre[s],*var1)
            reduce_chi2,aic,bic = check_module(n[s],modle_data,errors[s],i)
            differance = abs(reduce_chi2 - 1)
            if best_chi2 == -1 or best_chi2 > differance:
                best_chi2 = differance
                best_chi2_data = (var1,genlog_func,reduce_chi2,aic,bic)
            if best_aic == -1 or best_aic > aic:
                best_aic = aic
                best_aic_data = (var1,genlog_func,reduce_chi2,aic,bic)
            if best_bic == -1 or best_bic > bic:
                best_bic = bic
                best_bic_data = (var1,genlog_func,reduce_chi2,aic,bic)
        except Exception as e:
            print(e)
            print(i)
    for i in range(2,9,2):
        try:
            test = np.ones(i)
            var1 = curve_fit(genpow_func,bin_centre[s],n[s],sigma = errors[s],p0 = test)[0]
            modle_data = genpow_func(bin_centre[s],*var1)
            reduce_chi2,aic,bic = check_module(n[s],modle_data,errors[s],i)
            differance = abs(reduce_chi2 - 1)
            if best_chi2 == -1 or best_chi2 > differance:
                best_chi2 = differance
                best_chi2_data = (var1,genpow_func,reduce_chi2,aic,bic)
            if best_aic == -1 or best_aic > aic:
                best_aic = aic
                best_aic_data = (var1,genpow_func,reduce_chi2,aic,bic)
            if best_bic == -1 or best_bic > bic:
                best_bic = bic
                best_bic_data = (var1,genpow_func,reduce_chi2,aic,bic)
        except Exception as e:
            print(e)
            print(i)
    dictionary = {}
    dictionary["best_chi2_data"] = best_chi2_data
    dictionary["best_aic_data"] = best_aic_data
    dictionary["best_bic_data"] = best_bic_data
    best_results[pointer] = dictionary
def check_module(true_data,modle_data,errors,params):
        chi2 = get_chi2(true_data,modle_data,errors)
        reduce_chi2 = chi2 /(true_data.shape[0]-params)
        aic = chi2 + 2 * params
        bic = chi2 + params * np.log(true_data.shape[0])
        return reduce_chi2,aic,bic
def get_filedata(pointer):
    # get data from file
    with h5py.File(file_directory,"r") as file:
        #file.visit(general_attributes)
        #get_attributes("JWST/NIRCam.F356W/fiducial/N")
        x_data = file[f"JWST/NIRCam.{pointer}/fiducial/N"][:]
        y_data =  file[f"JWST/NIRCam.{pointer}/fiducial/bin_centres"][:]
        error_neg = file[f"JWST/NIRCam.{pointer}/fiducial/N_l"][:]
        error_pos = file[f"JWST/NIRCam.{pointer}/fiducial/N_u"][:]
    return x_data,y_data,error_neg,error_pos

for pointer in data_pointer:
    n,bin_centre,error_neg,error_pos = get_filedata(pointer)
    # curve fit doesnt take errors with differnant positvive/negtive values
    # so avarges are taken and infinate results are identified to be not included
    #errors = np.vstack([error_neg,error_pos])
    s=np.isfinite(n)
    errors = (error_pos + error_neg)/2

    find_func(bin_centre,s,errors,n,best_results,pointer)

chosen_data = []

print("here is the data")
for x in best_results.keys():
    for type in best_results[x].keys():
        print(f"best data for {type} is {best_results[x][type][1]} with reduce_chi2:{best_results[x][type][2]}  aic:{best_results[x][type][3]}  bic:{best_results[x][type][4]}")
    choose = input("input chi2,aic or bic to select best data: ")
    chosen_data.append(best_results[x][f"best_{choose}_data"])
# smooth the curve by using multiple x values
for data in chosen_data:
    data_type = data[1]
    best_data = data[0]
    print(genpol_print(*best_data))
    x = np.linspace(15.5,30,100)
    plt.plot(x,data_type(x,*best_data))

# plot the data
#plt.errorbar(bin_centre,n,yerr=errors)
plt.xlabel("magnitude")
plt.ylabel("galazys per unit area")
#print(genpol_print(*best_data))


plt.show()

