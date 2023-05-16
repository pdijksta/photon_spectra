import os
import subprocess
from PyQt5 import QtCore, QtWidgets

def send_to_desy_elog(author, title, severity, text, elog, image=None):
    """
    Send information to a supplied electronic logbook.
    Author: Christopher Behrens (DESY)
    """

    # The DOOCS elog expects an XML string in a particular format. This string
    # is beeing generated in the following as an initial list of strings.
    succeded = True  # indicator for a completely successful job
    # list beginning
    elogXMLStringList = ['<?xml version="1.0" encoding="ISO-8859-1"?>', '<entry>']

    # author information
    elogXMLStringList.append('<author>')
    elogXMLStringList.append(author)
    elogXMLStringList.append('</author>')
    # title information
    elogXMLStringList.append('<title>')
    elogXMLStringList.append(title)
    elogXMLStringList.append('</title>')
    # severity information
    elogXMLStringList.append('<severity>')
    elogXMLStringList.append(severity)
    elogXMLStringList.append('</severity>')
    # text information
    elogXMLStringList.append('<text>')
    elogXMLStringList.append(text)
    elogXMLStringList.append('</text>')
    # image information
    if image:
        try:
            # encodedImage = base64.b64encode(image)
            elogXMLStringList.append('<image>')
            elogXMLStringList.append(image)
            elogXMLStringList.append('</image>')
        except:  # make elog entry anyway, but return error (succeded = False)
            succeded = False
    # list end
    elogXMLStringList.append('</entry>')
    # join list to the final string
    elogXMLString = '\n'.join(elogXMLStringList)
    # open printer process
    try:
        lpr = subprocess.Popen(['/usr/bin/lp', '-o', 'raw', '-d', elog],
                               stdin=subprocess.PIPE, stdout=subprocess.PIPE)
        # send printer job
        lpr.communicate(elogXMLString.encode('utf-8'))
    except:
        succeded = False
    return succeded

def logbook(widget, text=""):
    screenshot = get_screenshot(widget)
    res = send_to_desy_elog(author='Dr. Snail', title='OrbitSnailScan', severity='INFO', text=text, elog='xfellog', image=screenshot)
    if not res:
        print('error during eLogBook sending')

def log_screen(save_func, widget, auto_comment=""):
    filename, comment, ok = save_func()
    text = 'Data is saved in %s' % os.path.abspath(filename)
    if ok:
        text = comment + "\n" +"\n" + text
    text = auto_comment + text
    logbook(widget, text=text)

def get_screenshot(window_widget):
    screenshot_tmp = QtCore.QByteArray()
    screeshot_buffer = QtCore.QBuffer(screenshot_tmp)
    screeshot_buffer.open(QtCore.QIODevice.WriteOnly)
    widget = QtWidgets.QWidget.grab(window_widget)
    widget.save(screeshot_buffer, "png")
    return screenshot_tmp.toBase64().data().decode()

def dialog(parent):
    dlg = QtWidgets.QInputDialog(parent)
    dlg.setInputMode(QtWidgets.QInputDialog.TextInput)
    dlg.setLabelText("Comment :")
    dlg.resize(400, 100)
    ok = dlg.exec_()
    comment = dlg.textValue()
    return ok, comment

