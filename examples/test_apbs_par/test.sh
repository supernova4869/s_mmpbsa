rm cost*
for i in 1 4 8 16 32 64 128
do
    export OMP_NUM_THREADS=$i
    for j in `seq 1 3`
    do
        echo "Running with $i kernels the "$j"th time"
        (time apbs test.apbs) 1>>/dev/null 2>> cost_$i.txt
    done
done

echo -e "cores,Time (s)" > cost.csv
for i in 1 4 8 16 32 64 128
do
    for v in $(awk '/real/ {print $2;}' cost_$i.txt | awk -F "m" '{print $1 * 60 + $2;}')
    do
        echo -e "$i,$v" >> cost.csv
    done
done