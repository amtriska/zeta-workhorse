from zeta_workhorse import load_gml, tropical_zeta

try:
    data = load_gml("karate.gml", data_type = float)
    tropical_result = tropical_zeta(data, 5)
    print("The sequence of minimum cycle costs up to length 5 is: ", tropical_result)
    input("Press Enter to continue.\n")
except Exception as error:
    print("Calculation not execute successfully. See error message below.")
    print(error)
    input("Press Enter to continue.\n")