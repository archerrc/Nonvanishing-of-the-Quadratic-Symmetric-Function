︠d718b6a7-46e3-405e-8022-da45caecca8bs︠
def P2k(k,t,m): #here k is half of the weight
    Ptot=0
    for j in range (0,k):
        #this is the combinatorial formula given in Chiriac and Jorza's paper "The trace of T2 takes no repeated values"
        Ptot= (-1)^j*binomial(2*k-2-j,j)*m^j*t^(2*k-2-2*j)+Ptot 
    return(Ptot)

#for TrT2 and TrT4, these are from the version of the Eichler-Selberg trace formula found in Zagier's appendix to Lang's "Intorduction to modular forms"
def TrT2(k): #k is half of the weight
    Trace2=-1/2*P2k(k,0,2)-P2k(k,1,2)-1/2*P2k(k,2,2)-1
    return(Trace2)

def TrT4(k):
    Trace4=-3/4*P2k(k,0,4)-2*P2k(k,1,4)-4/3*P2k(k,2,4)+1/12*P2k(k,4,4)-2^(2*k-2)-1
    return(Trace4)

def EigenSum(k): #This is sum(a(i)a(j))
    A=TrT2(k)^2
    B=TrT4(k)
    dk=floor(k/6)
    if k%6==1:
        dk=dk-1
    Eigens=1/2*(A-B)-2^(2*k-2)*dk
    return(Eigens)
︡8757e7b8-5cd4-4775-ac8a-693432015d23︡{"done":true}
︠13f08dd1-033a-4a74-aeca-b00ff18ca53as︠
for i in range (1,68):
    if EigenSum(i)==0:
            print('zero at ', i)
    if i==67:
            print('Done')
︡2ce85127-1e4e-4d04-9501-0756840ffbd0︡{"stdout":"Done"}︡{"stdout":"\n"}︡{"done":true}
︠11efe5f0-31aa-4a71-9efb-5f6b0ebae8ab︠









