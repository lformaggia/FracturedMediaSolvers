#!/bin/bash

DATE=`date +%x`
DATE_TAG=`date +%Y%m%d`
TAG_NAME="Release_Kafres-Upscaling_${DATE_TAG}"
RELEASE_FILE="RELEASE"

echo ${TAG_NAME} > ${RELEASE_FILE}
git commit -m "Updated RELEASE file for version ${TAG_NAME}" ${RELEASE_FILE}
echo "Updated RELEASE file and commited the change"
git tag -a "${TAG_NAME}" -m "Version released to ENI on ${DATE}"
echo "New tag ${TAG_NAME} created."
