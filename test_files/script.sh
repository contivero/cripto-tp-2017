../bin/bmpsss -d --secret Albertssd.bmp -w 300 -h 300 -k 8 --dir g7
mv shadow* shades
../bin/bmpsss -r --secret o.bmp -k 8 -w 300 -h 300 --dir shades

../bin/bmpsss -r --secret a.bmp -k 8 -w 450 -h 300 --dir g7
