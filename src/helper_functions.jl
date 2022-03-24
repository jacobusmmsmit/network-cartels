norme(x) = √((x[1])^2 + (x[2])^2)
diste(x, y) = norme(x .- y)
disth(x, y) = acosh(1 + (2 * diste(x, y)^2) / ((1 - norme(x)^2) * (1 - norme(y)^2)))

# UnitRangeWithOut
ur_wo(ur, without) = Iterators.flatten((ur.start:(without-1), (without+1):ur.stop))