from selenium.webdriver import Chrome, ChromeOptions
from selenium.webdriver.common.by import By
from selenium.webdriver.support.ui import WebDriverWait
from selenium.webdriver.support import expected_conditions as EC

mail = "kmatsumoto@gen-info.osaka-u.ac.jp"
pwd  = "AFCryfx5"
driver_path = '/Users/kaoru/chromedriver/chromedriver'

# ドライバを得る
options = ChromeOptions()
#options.add_argument('--headless')
browser = Chrome(executable_path=driver_path ,options=options)

# ログインページにアクセス
url_login = "https://app.asana.com/-/login"
browser.get(url_login)
print("ログインページにアクセスしました")

# テキストボックスに文字を入力 --- (※3)
e = browser.find_element_by_name("e")
e.clear()
e.send_keys(mail)
e = browser.find_element_by_name("p")
e.clear()
e.send_keys(pwd)

#ログインボタンをクリック
login_btn = browser.find_element_by_class_name('LoginEmailPasswordForm-logInButton')
login_btn.click()
print("情報を入力してログインボタンを押しました")

# ページのロード完了まで待機
WebDriverWait(browser, 15).until(
    EC.presence_of_element_located((By.XPATH, "//*[contains(text(), 'Asana へようこそ。')]"))
)

# マイタスクのURLを得て移動
mytask     = browser.find_element_by_class_name('SidebarTopNavLinks-myTasksButton')
mytask_url = mytask.get_attribute('href')
browser.get(mytask_url)


task_cells = browser.find_elements_by_class_name('SpreadsheetGridTaskNameCell')
task_names = [i.get_attribute('aria-label') for i in task_cells]

browser.close()
