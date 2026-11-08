R<z> := PolynomialRing(Integers());

print "Form 35629";
print Solutions(Thue(z^3 + 18 * z^2 - 36 * z - 72), 1 : Exact := true);

print "Form 50870";
print Solutions(Thue(z^3 + 27 * z^2 - 81 * z - 243), 1 : Exact := true);

print "Form 78530";
print Solutions(Thue(z^3 + 45 * z^2 - 69 * z - 97), 1 : Exact := true);

print "Form 85192";
print Solutions(Thue(z^3 + 48 * z^2 + 24 * z - 120), 1 : Exact := true);

print "Form 217159";
print Solutions(Thue(3 * z^3 + 30 * z^2 - 84 * z - 32), 1 : Exact := true);

print "Form 223187";
print Solutions(Thue(3 * z^3 + 39 * z^2 - 15 * z - 83), 1 : Exact := true);

print "Form 407524";
print Solutions(Thue(7 * z^3 + 57 * z^2 - 3 * z - 29), 1 : Exact := true);

print "Form 441326";
print Solutions(Thue(8 * z^3 + 48 * z^2 - 66 * z - 17), 1 : Exact := true);

print "Form 447517";
print Solutions(Thue(8 * z^3 + 60 * z^2 - 12 * z - 37), 1 : Exact := true);

print "Form 460471";
print Solutions(Thue(9 * z^3 + 18 * z^2 - 36 * z - 8), 1 : Exact := true);

print "Form 594271";
print Solutions(Thue(15 * z^3 + 12 * z^2 - 96 * z + 8), 1 : Exact := true);

print "All eleven equations solved.";
