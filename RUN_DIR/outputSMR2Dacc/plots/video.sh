#!/bin/sh

#while getopts i:o: flag
#do
#        case "${flag}" in
#                i) input=${OPTARG};;
#                o) output=${OPTARG};;
#        esac
#done

#echo "$input"
echo $1
#ffmpeg -framerate 5 -i cbd.out1.0%04d_rho.png -vcodec libx264 -crf 18 rho.mp4
#ffmpeg -framerate 5 $input -vcodec libx264 -crf 18 $output
ffmpeg -framerate 5 $1 -vcodec libx264 -crf 18 $2
