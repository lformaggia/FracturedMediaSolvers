#! /bin/bash

DATE=`date +%x`
DATE_TAG=`date +%Y%m%d`
TAG_NAME="Release_Kafres-Upscaling_${DATE_TAG}"
RELEASE_FILE="RELEASE"

git tag -a "${TAG_NAME}" -m "Version released to ENI on ${DATE}"
echo $TAG_NAME > $RELEASE_FILE
echo -e "New tag ${TAG_NAME} created.\n\
Remember to push the tag and the ${RELEASE_FILE} file.\n\
Type the following code to push the tag:\n\
git push origin ${TAG_NAME}"
