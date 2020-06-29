# Mendelian Inheritance

def autosomal(mal, fem, exp):
    '''
    Autosomal dominant ("D"):
    - each affected child has an affected parent
    - occurs in every generation
    Autosomal recessive ("R"):
    - both parents of affected child are carriers
    - doesn't usually occur in every generation
    '''
    homod = mal.upper() # homozygous dominant 
    homor = mal.lower() # homozygous recessive
    hetero = homod[0] + homor[0] # heterozygous

    # offspring genotype percentages
    # pct = [AA, Aa, aa]
    if mal == homod and fem == homod:
        pct = [1, 0, 0]
    elif mal == hetero and fem == hetero: 
        pct = [0.25, 0.5, 0.25]
    elif mal == homor and fem == homor:
        pct = [0, 0, 1]
    elif mal == homod and fem == homor or mal == homor and fem == homod:
        pct = [0, 1, 0]
    elif mal == homod and fem == hetero or mal == hetero and fem == homod:
        pct = [0.5, 0.5, 0]
    elif mal == hetero and fem == homor or mal == homor and fem == hetero:
        pct = [0, 0.5, 0.5]

    # phenotype expression percentages 
    if exp.upper() == "D":
        y = pct[0] + pct[1]
        n = pct[2]
    elif exp.upper() == "R":
        y = pct[1] + pct[2]
        n = pct[0]  
    
    print((homod + ": {}\n" + hetero + ": {}\n" + homor + ": {}\n").format(*pct))
    print("Express:  %.2f\nNo Express: %.2f" % (y,n))

def x_linked(mal, fem, exp):
    '''
    X-linked Dominant ("D"):
    - females more frequently affected
    - can have both affected males and females in the same generation
    X-linked Recessive ("R"):
    - males more frequently affected
    - affected males often present in each generation
    '''
    homod = "XX" # homozygous dominant female
    homor = "xx" # homozygous recessive female
    hetero = "Xx" # heterozygous female
    maled = "XY" # dominant male
    maler = "xY" # recessive male

    # offspring genotype percentages
    # pct = [XX, Xx, xx, XY, xY]
    if mal == maled:
        if fem == homod:
            pct = [0.5, 0, 0, 0.5, 0]
        elif fem == hetero:
            pct = [0.25, 0.25, 0, 0.25, 0.25]
        elif fem == homor:
            pct = [0, 0.5, 0, 0, 0.5]
    elif mal == maler:
        if fem == homod:
            pct = [0, 0.5, 0, 0.5, 0]
        elif fem == hetero:
            pct = [0, 0.25, 0.25, 0.25, 0.25]
        elif fem == homor:
            pct = [0, 0, 0.5, 0, 0.5]
    
    # phenotype expression percentages
    if exp.upper() == "D":
        y = pct[0] + pct[1] + pct[3]
        n = pct[2] + pct[4]
    elif exp.upper() == "R":
        y = pct[1] + pct[2] + pct[4]
        n = pct[0] + pct[3]

    print((homod + ": {}\n" + hetero + ": {}\n" + homor + ": {}\n" + maled + ": {}\n" + maler + ": {}\n").format(*pct))
    print("Express:  %.2f\nNo Express: %.2f" % (y,n))

def y_linked(mal, fem, exp):
    '''
    Y-linked Dominant ("D") and Y-linked Recessive ("R") can only be passed from father to son
    '''
    homo = "XX" # female
    maled = "XY" # dominant male
    maler = "Xy" # recessive male
        
    # offspring genotype percentages
    # pct = [XX, XY, Xy]
    if mal == maled:
        pct = [0.5, 0.5, 0]
    elif mal == maler:
        pct = [0.5, 0, 0.5]

    # phenotype expression percentages
    if exp.upper() == "D":
        y = pct[1]
        n = pct[0] + pct[2]
    elif exp.upper() == "R":
        y = pct[2]
        n = pct[0] + pct[1]

    print((homo + ": {}\n" + maled + ": {}\n" + maler + ": {}\n").format(*pct))
    print("Express:  %.2f\nNo Express: %.2f" % (y,n))








    
