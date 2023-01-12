const fs = require('fs');
const xml2js = require('xml2js');
const path = require('path');
const glob = require("glob")

const data_path = path.resolve('./download');

var files = glob.sync(data_path + '/*.rm5');

let j = 4;

var parser = new xml2js.Parser();
fs.readFile(files[j], function(err, data) {
    parser.parseString(data, function (err, result) {
        console.dir(result);
        let authors = result['COCHRANE_REVIEW']['COVER_SHEET'][0]['CREATORS'][0]['PERSON'].map(x => x['FIRST_NAME'] + ' ' + x['MIDDLE_INITIALS'] + ' ' + x['LAST_NAME']);
        let title = result['COCHRANE_REVIEW']['COVER_SHEET'][0]['TITLE'][0].replace(/\[.*/i,'').trim();
        console.log(title);
        console.log(files[j]);
        console.log('');
    });
});