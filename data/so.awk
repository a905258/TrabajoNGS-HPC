NR==1{print $0, "   num.mark.\n"}
NR>2{printf("%s\t%6.0f\n", $0, $3-$2)}

