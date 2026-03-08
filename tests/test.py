from zeta_workhorse import load_gml, load_csv, tropical_trace, ihara_zeta_poly

try:
    print("Calculating...")    
    trop_data = load_gml("karate_weighted.gml", data_type = float, weight_field = "value")
    ihara_data = load_csv("barbell.csv", data_type = float)
    tropical_result = tropical_trace(trop_data, 10, min)
    ihara_result = ihara_zeta_poly(ihara_data)
    print("The sequence of minimum cycle costs up to length 10 for the weighted Karate Club dataset is: ")
    print(tropical_result)
    print("The Ihara zeta polynomial of the barbell graph dataset is: ")
    print(ihara_result)
    input("Press Enter to continue.\n")
except Exception as error:
    print("Calculation not executed successfully. See error message below.")
    print(error)
    input("Press Enter to continue.\n")