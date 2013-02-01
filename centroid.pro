FUNCTION centroid, input

size = size(input, /Dimensions)

totalMass = total(input)

xcom = total(total(input,2) * IndGen(size[0])) / totalMass
ycom = total(total(input,1) * IndGen(size[1])) / totalMass

RETURN, [xcom,ycom]
END
