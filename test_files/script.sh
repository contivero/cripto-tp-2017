../bin/bmpsss -d --secret Albertssd.bmp -w 300 -h 300 -k 8 --dir imgs_450x300
mkdir -p shades-tmp
mv shadow*.bmp shades-tmp
mkdir -p outputs
../bin/bmpsss -r --secret outputs/distributed_and_recovered.bmp -k 8 -w 300 -h 300 --dir shades-tmp

../bin/bmpsss -r --secret outputs/output1.bmp -k 8 -w 450 -h 300 --dir imgs_450x300
../bin/bmpsss -r --secret outputs/output2.bmp -k 8 -w 300 -h 450 --dir imgs_300x450
../bin/bmpsss -r --secret outputs/output3.bmp -k 8 -w 300 -h 300 --dir imgs_300x300
../bin/bmpsss -r --secret outputs/output4.bmp -k 2 -w 600 -h 398 --dir imgs_600x1593
