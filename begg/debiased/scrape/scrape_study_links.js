const puppeteer = require('puppeteer');
const fs = require('fs');
const umass_login = require('/Users/haben/OneDrive - University of Massachusetts/scripts/umass_login');
const save_file = './study_links.txt';


(async () => {
    const browser = await puppeteer.launch({
        headless: false,
        args: ['--window-size=1920,1040'],
        slowMo: 100
    });
    const page = await browser.newPage();
    await page.setViewport({ width: 1920, height: 1040 });
    await page.goto('http://silk.library.umass.edu/login?url=https://www.cochranelibrary.com/advanced-search');

    await umass_login(page);

    await page.waitForSelector('#searchText0');
    await page.$eval('#searchText0', el => el.value = 'meta-analysis');
    await page.waitForSelector('button[type=submit]');
    await page.click('button[type=submit]');


    await page.waitForSelector('#orderBy');
    await page.select('#orderBy', 'displayDate-true');

    while (true) {
        // let page_links = await page.evaluate(function () {
        //     let links = [];
        //     document.querySelectorAll(".result-title").forEach(function (e, i, j) {
        //         links.push(e.querySelector('a').getAttribute('href'))
        //     });

        //     return links;
        // });
        let page_links = await page.$$eval(
            '.result-title > a',
            as => as.map(a => a.href)
        );

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
