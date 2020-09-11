BEGIN{
    n=0; fi="";
}

{
    if ($0 ~ /\*FG\*/) {
        if (n>=5) {
            print s;
            print fi;
        }
        s=$0;
        n=0;
        fi="";
    } else {
        if ($NF<0) {
            n=n+1;
            fi= sprintf("%s%s\n", fi, $0);
        }
    }
}
