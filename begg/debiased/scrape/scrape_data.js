const puppeteer = require('puppeteer');
const fs = require('fs');
//const path = require('path');
const https = require('https');
const umass_login = require('/Users/haben/OneDrive - University of Massachusetts/scripts/umass_login');
const path = require('path');
const data_links_file = 'data_links.txt';
const umass_cochran_link = 'https://www-cochranelibrary-com.silk.library.umass.edu/';

let links = fs.readFileSync(data_links_file, 'utf8');
links = links.split('\n').map(s => s.split('###')[1]);
links = [...new Set(links)];


(async () => {
    const browser = await puppeteer.launch({
        headless: false,
        args: ['--window-size=1920,1040'],
        slowMo: 100
    });
    const page = await browser.newPage();
    await page.setViewport({ width: 1920, height: 1040 });
    await page.goto(umass_cochran_link);
    umass_login(page);
    await page.waitForTimeout(1000);
    console.log('dsljf');

    for (let j = 0; j < links.length; j++) {
        console.log(links[j]);
        await page._client.send('Page.setDownloadBehavior', {behavior: 'allow', downloadPath: path.resolve('./download')});
        try {
            await page.waitForTimeout(3000);
            await page.goto(links[j]);

        // await page.goto(umass_cochran_link);

        } catch(error) {
            
        }
        // const file = fs.createWriteStream('./download/test.txt');
        // const request = https.get('https://www-cochranelibrary-com.silk.library.umass.edu/cdsr/doi/10.1002/14651858.CD011841.pub2/media/CDSR/CD011841/table_n/CD011841StatsDataOnly.rm5?content-disposition=attachment&mime-type=application/octet-stream', function (response) {
        //     response.pipe(file);
        // });
    }
})();
