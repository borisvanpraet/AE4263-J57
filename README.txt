	Git Tutorial
Alle Git taken moeten worden uitgevoerd vanuit de command line of de Git Bash/GUI.
Nadat je de repository hebt gecloned, moet je altijd eerst git init gebruiken voordat je andere taken volbrengt.

- Clone repository:
git clone https://github.com/Floris89/AE4263-J57

- File toevoegen aan repository
Om alle files in een folder toe te voegen: git add .
Om 1 file toe te voegen: git add "file"  (file vervangen voor de file die je wil toevoegen)

- Commits
git commit -a -m "Message" (Message vervangen voor een commit bericht)

- Branching
Om een nieuwe branch te maken: git branch my-new-branch	(my-new-branch vervangen voor een goede naam)
Om alle branches te zien: git branch
Om naar een branch te gaan: git checkout my-new-branch (my-new-branch vervangen voor de branch naam die je gebruikt)

- Push
Als je voor het eerst iets pusht: git push -u origin main
Daarna: git push origin my-branch (my-branch vervangen voor de branch waarnaar je pusht)
Om een nieuwe branch te pushen: git push -u origin my-new-branch (my-new-branch vervangen voor de branch name die je gebruikt)

- Pull
git fetch origin