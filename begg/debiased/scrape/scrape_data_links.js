const puppeteer = require('puppeteer');
const fs = require('fs');
const https = require('https');
const umass_login = require('/Users/haben/OneDrive - University of Massachusetts/scripts/umass_login');
const study_links_file = 'study_links.txt';
const save_file = './data_links.txt';
const umass_cochran_link = 'https://www-cochranelibrary-com.silk.library.umass.edu/';
const base_path = 'https://www-cochranelibrary-com.silk.library.umass.edu';

let links = fs.readFileSync(study_links_file, 'utf8');
links = links.split('\n');
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
    await umass_login(page);

    for (let j = 0; j < links.length; j++) {
        try {
            console.log(j);
            await page.waitForTimeout(1000);
            await page.goto(links[j]);
            await page.waitForSelector('.download-stats-data-link', {timeout: 5000});
            let data_links = await page.$$eval(
                '.download-stats-data-link > a',
                as => as.map(a => a.href)
            );

            out_string = links[j] + '###' + data_links[0]  + '\n';
            fs.appendFileSync(save_file, out_string);
        } catch (error) {
            console.log(error);
        }
    };



    // await page.waitForSelector('.cdsr-nav-link download-stats-data-link');
    // '.cdsr-nav-link download-stats-data-link > a'
    // const file = fs.createWriteStream("tt.txt");
    // const request = http.get("http://i3.ytimg.com/vi/J---aiyznGQ/mqdefault.jpg", function(response) {
    //   response.pipe(file);

    await page.$eval('#searchText0', el => el.value = 'meta-analysis');
    await page.waitForSelector('button[type=submit]');
    await page.click('button[type=submit]');

    // await page.evaluate(_ => {
    //     window.scrollBy(0, window.innerHeight);
    //     //window.scrollBy({top:300});
    //   });   
    //   await page.waitForTimeout(20000);

    await page.waitForSelector('#orderBy');
    await page.select('#orderBy', 'displayDate-true');

    while (true) {
        let page_links = await page.evaluate(function () {
            let links = [];
            document.querySelectorAll(".result-title").forEach(function (e, i, j) {
                links.push(e.querySelector('a').getAttribute('href'))
            });
            return links;
        });

        console.log(page_links);
        await page.waitForSelector('.pagination-next-link');
        await page.waitForTimeout(3000);
        await page.click(".pagination-next-link");
        await page.waitForSelector('.result-title');
        await page.waitForTimeout(3000);

        out_string = page_links.join('\n').concat('\n');
        fs.appendFileSync(save_file, out_string);
        // await page.waitForTimeout(1000);
    };

}
)();
