for i in {1..6}
do
  sleep $((2 * $i)) &&
  ../bw-s1p2-cpu > 36-$i.out &
done

