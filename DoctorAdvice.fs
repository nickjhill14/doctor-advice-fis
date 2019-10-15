# A fuzzy inference system (FIS) for advising a doctor whether a patient should be referred to a hospital for emergency investigations
# Two biomedical inputs: The patient’s temperature and the severity of headache

# Creating a new FIS
fis = newfis("Referrels", defuzzMethod = "bisector")
#fis = newfis("Referrels", defuzzMethod = "som")

# INPUTS
# Temperature is represented by degrees celcius       
fis = addvar(fis, "input", "Temperature", c(34, 39))
# Headache severity represented by an arbitrary scale from 0 to 10
fis = addvar(fis, "input", "Headache", c(0, 10))

# OUTPUT
# Urgency (or severity) of the situation is represented by an arbitrary scale from 0 to 100
fis = addvar(fis, "output", "Urgency", c(0, 100))

# MEMBERSHIP FUNCTIONS
# Temperature linguistic variables
fis = addmf(fis, "input", 1, "Hypothermic", "trimf", c(-1e100, 34, 35))
fis = addmf(fis, "input", 1, "Low", "trimf", c(34, 35.5, 36.5))
fis = addmf(fis, "input", 1, "Normal", "trimf", c(35.5, 36.5, 37.5))
fis = addmf(fis, "input", 1, "High", "trimf", c(36.5, 37.5, 39))
fis = addmf(fis, "input", 1, "Hyperthermic", "trimf", c(38, 39, +1e100))

# Headache severity linguistic variables
fis = addmf(fis, "input", 2, "No Pain", "trimf", c(-1e100, 0, 2.5))
fis = addmf(fis, "input", 2, "Mild", "trimf", c(0, 2.5, 5))
fis = addmf(fis, "input", 2, "Moderate", "trimf", c(2.5, 5, 7.5))
fis = addmf(fis, "input", 2, "Severe", "trimf", c(5, 7.5, 10))
fis = addmf(fis, "input", 2, "Very Severe", "trimf", c(7.5, 10, +1e100))

# Urgency (or severity) of the situation linguistic variables
fis = addmf(fis, "output", 1, "Non-Urgent", "gaussmf", c(10, 0))
fis = addmf(fis, "output", 1, "Less Urgent", "gaussmf", c(10, 25))
fis = addmf(fis, "output", 1, "Urgent", "gaussmf", c(10, 50))
fis = addmf(fis, "output", 1, "Emergent", "gaussmf", c(10, 75))
fis = addmf(fis, "output", 1, "Resuscitation", "gaussmf", c(10, 100))

# PLOTTING MEMBERSHIP FUNCTIONS
#plotmf(fis, "input", 1, NULL, 0, "Temperature (Celcius)", "Degree of Membership", 'Input Membership Functions for the Linguistic Variable "Temperature"')
#plotmf(fis, "input", 2, NULL, 0, "Headache Severity (Ranked from 0 to 10)",  "Degree of Membership", 'Input Membership Functions for the Linguistic Variable "Headache"')
#plotmf(fis, "output", 1, NULL, 0, "Urgency", "Degree of Membership", 'Output Membership Functions for the Linguistic Variable "Urgency"')

# OPERATOR
fuzzy.t(min, fis)

# RULES
ruleList = rbind(c(1, 1, 5, 1, 1), c(1, 2, 5, 1, 1), c(1, 3, 5, 1, 1), c(1, 4, 5, 1, 1), c(1, 5, 5, 1, 1), c(2, 1, 2, 1, 1), c(2, 2, 3, 1, 1), c(2, 3, 3, 1, 1), c(2, 4, 4, 1, 1), c(2, 5, 5, 1, 1), c(3, 1, 1, 1, 1), c(3, 2, 1, 1, 1), c(3, 3, 2, 1, 1), c(3, 4, 3, 1, 1), c(3, 5, 5, 1, 1), c(4, 1, 2, 1, 1), c(4, 2, 3, 1, 1), c(4, 3, 3, 1, 1), c(4, 4, 4, 1, 1), c(4, 5, 5, 1, 1), c(5, 1, 5, 1, 1), c(5, 2, 5, 1, 1), c(5, 3, 5, 1, 1), c(5, 4, 5, 1, 1), c(5, 5, 5, 1, 1))
fis = addrule(fis, ruleList)

# TESTING
# Rule list testing
inputs = rbind(c(34, 0), c(34, 2.5), c(34, 5), c(34, 7.5), c(34, 10), c(35, 0), c(35.25, 2.5), c(35.25, 5), c(35.25, 7.5), c(35.25, 10), c(36.5, 0), c(36.5, .25), c(36.5, 5), c(36.5, 7.5), c(36.5, 10), c(37.75, 0), c(37.75, 2.5), c(37.75, 5), c(37.75, 7.5), c(37.75, 10), c(39, 0), c(39, 2.5), c(39, 5), c(39, 7.5), c(39, 10))
# Performance analysis testing
#inputs = rbind(c(36.5, 0), c(39, 10), c(34, 10), c(36.5, 10), c(39, 0), c(34, 0), c(36.5, 5), c(35, 5), c(38, 5), c(36.5, 2.5), c(36.5, 7.5))
evalfis(inputs, fis)
outputs = evalfis(inputs, fis)

# Prints out the FIS
showfis(fis)
# Prints out the rules
showrule(fis)

# Visualise the output surface plot for the two-input one-output FIS 
gensurf(fis)