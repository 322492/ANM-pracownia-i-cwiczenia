# autor: Kamil Tasarz

using Printf

#przygotowanie

two_to = (2.0 .^ (1:62)) #tablica kolejnych potęg dwójki
four_to = (4.0 .^ (1:31)) #tablica kolejnych potęg czwórki

function trapez(f, a, b, n, h)
    res = 0.0
    for i in 0:n
        if i==0 || i==n
            res += f(a + i*h) / 2.0
        else
            res += f(a + i*h)
        end
    end
    return res*h
end

function summa(to, f, a, h) # aktulanie nieużywana, można użyć zamiast trapez
    s = 0.0
    for i in 1:to
        s = s + f(a + (2.0 * i - 1.0) * h)
    end
    return s
end

# R_(n, m), gdzie n - wiersz, m - kolumna
# przy czym Julia numeruje tablice od 1 a nie od zera, więc wyraz R_(i, j) jest w miejscu M[i+1][j+1]

function Romberg_KT(a, b, f, ϵ, drukowac=false, exact=0.0)
    # początek i koniec przedzialu całki, funkcja, oczekiwany błąd, czy wypisać tablicę błędów oraz
    # dokładna wartość całki (potrzebna do wypisania tablicy błędów; jeśli nie podamy otrzymamy po prostu tablicę Romberga)
    # zwracana wynik to liczba wierszy i najlepsze uzyskane przybliżenie całki
    M = [] #tablica Romberga
    h = (b - a)
    R00 = h * (f(b)-f(a)) / 2.0
    push!(M, R00)
    i = 1 #numer wiersza
    h = (b - a)
    
    while true
        w = [] #aktualny wiersz
        h = h / 2.0
        # R = M[i][1]/2.0 + h*summa(1, two_to[i]/2, f, a, h) #nieco szybsze
        R = trapez(f, a, b, two_to[i], h)
        push!(w, R)
        
        for j in 1:1:i
           # R = w[j] + (w[j] - M[i][j])/(four_to[j] - 1.0)
           R = (four_to[j] * w[j] - M[i][j]) / (four_to[j] - 1.0) # drugi równoważny sposób
            push!(w, R)
        end
        push!(M, w)
        
        if(abs(M[i+1][1] - M[i][1]) < ϵ*abs(M[i+1][1])) #warunek z treści zadania: |R_(K, 0) - R_(K-1, 0)| < ϵ*|R_(K, 0)|
            break
        end
        
        i+=1
    end
            
    if drukowac
        #drukuj_ostatnie_wyrazy_tablicy(M)
        drukuj_tablice_bledow(M, exact)
    end
    
    return i, M[i+1][i+1]
end

     

function drukuj_tablice_bledow(M, exact)
    if exact == 0.0 # drukowanie po prostu tablicy Romberga
        for i in 1:1:length(M)
            @printf("%d: ", i-1)
            for j in 1:1:i
                @printf("%.5f ", M[i][j])
            end
            @printf("\n")
        end
    else #drukowanie błędów
        for i in 1:1:length(M)
            @printf("%d: ", i-1)
            for j in 1:1:i
                @printf("%.3e ", abs(exact - M[i][j]) / abs(exact))
            end
            @printf("\n")
        end
    end  
end

function drukuj_ostatnie_wyrazy_tablicy(M) # zamiast wypisywać całe wiersze można wypisać najlepsze przybliżenie z każdego
    for i in 1:1:length(M)
            @printf("%d: %f \n", i-1, M[i][i])
        end
end

#test1 - wielomian
a = -1.0
b = 1.0
f(x) = 2.0 * x^3 - 6.0/5.0 * x^2 + 0.17
exact = -0.46

print("epsilon i wynik funckji Romberg", "\n")
for ϵ in [1e-1, 1e-2, 1e-3, 1e-4, 1e-5, 1e-6, 1e-7, 1e-8]
    @printf("%.1e  ", ϵ)
    print(Romberg_KT(a, b, f, ϵ, false, exact), "\n")
end

#test2 - funkcja wykładnicza
a = -1.0
b = 1.0
f(x) = ℯ^x
exact = 2.350402387287602913

print("epsilon i wynik funckji Romberg", "\n")
for ϵ in [1e-1, 1e-2, 1e-3, 1e-4, 1e-5, 1e-6, 1e-7, 1e-8]
    @printf("%.1e  ", ϵ)
    print(Romberg_KT(a, b, f, ϵ, false, exact), "\n")
end

#test3 - funkcja Rungego
a = -1.0
b = 1.0
f(x) = 1.0 / (1.0 + 25.0 * x^2)
exact = 0.549360306778006344344508770577

print("epsilon i wynik funckji Romberg", "\n")
for ϵ in [1e-1, 1e-2, 1e-3, 1e-4, 1e-5, 1e-6, 1e-7, 1e-8]
    @printf("%.1e  ", ϵ)
    print(Romberg_KT(a, b, f, ϵ, false, exact), "\n")
end

#test4 - z treści zadania
a = -1.0
b = 1.0
f(x) = 1.0 / (x^4 + x^2 + 0.9)
exact = 1.58223

print("epsilon i wynik funckji Romberg", "\n")
for ϵ in [1e-1, 1e-2, 1e-3, 1e-4, 1e-5, 1e-6, 1e-7, 1e-8]
    @printf("%.1e  ", ϵ)
    print(Romberg_KT(a, b, f, ϵ, false, exact), "\n")
end

#test5 - z treści zadania
a = 0.0
b = 1.0
f(x) = 1.0 / (1.0 + x^4)
exact = 0.866972987339911037573995163882870713652175367

print("epsilon i wynik funckji Romberg", "\n")
for ϵ in [1e-1, 1e-2, 1e-3, 1e-4, 1e-5, 1e-6, 1e-7, 1e-8]
    @printf("%.1e  ", ϵ)
    print(Romberg_KT(a, b, f, ϵ, false, exact), "\n")
end

#test6 - z treści zadania
a = 0.0
b = 1.0
f(x) = 2.0 / (2.0 + sin(10*pi*x))
ϵ = 1e-8
exact = 1.1547005383792515290182975610039149 # 2/sqrt(3)

print("epsilon i wynik funckji Romberg", "\n")
for ϵ in [1e-1, 1e-2, 1e-3, 1e-4, 1e-5, 1e-6, 1e-7, 1e-8]
    @printf("%.1e  ", ϵ)
    print(Romberg_KT(a, b, f, ϵ, false, exact), "\n")
end

#test drukowania tablicy błędów
a = -1.0
b = 1.0
ϵ = 1e-6
f(x) = sqrt(sqrt(sqrt(1.0 - x))) # (1-x)^(1/8)
exact = 1.9386804136271247274791300546857030

Romberg_KT(a, b, f, ϵ, true, exact)
